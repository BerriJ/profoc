// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <misc.h>
#include <loss.h>
#include <splines2.h>
#include <splines.h>

#include <RcppArmadillo.h>
#include <progress.hpp>

using namespace arma;

void online_learning_core_mv(
    const unsigned int &T,
    const unsigned int &D,
    const unsigned int &P,
    const unsigned int &K,
    const unsigned int &T_E_Y,
    const unsigned int &start,
    const unsigned int &lead_time,
    const mat &y,
    const field<cube> &experts,
    const vec &tau,
    const std::string loss_function,
    const std::string method,
    const double &loss_parameter,
    const bool &loss_gradient,
    field<cube> &weights,
    field<cube> &weights_tmp,
    field<cube> &predictions_tmp,
    cube &predictions,
    field<cube> &past_performance,
    vec &opt_index,
    const mat &param_grid,
    mat &chosen_params,
    field<cube> &R,
    field<cube> &R_reg,
    field<cube> &V,
    field<cube> &E,
    field<cube> &eta,
    field<mat> &hat_mats,
    field<sp_mat> &basis_mats,
    field<sp_mat> &basis_mats_mv,
    field<cube> &beta0field,
    field<cube> &beta,
    mat &cum_performance,
    const double &forget_past_performance,
    mat &tmp_performance,
    const bool &allow_quantile_crossing,
    cube &loss_for,
    field<cube> &loss_exp,
    const cube &loss_array,
    const cube &regret_array,
    Progress &prog)
{
  for (unsigned int t = start; t < T; t++)
  {
    weights(t).set_size(D, P, K);
    past_performance(t).set_size(D, P, param_grid.n_rows);

    for (unsigned int d = 0; d < D; d++)
    {
      // Save final weights weights_tmp
      weights(t).row(d) = weights_tmp(opt_index(t)).row(d);

      // Store expert predictions temporarily
      mat experts_mat = experts(t).row(d);

      // Predictions using different parameter values
      for (unsigned int x = 0; x < param_grid.n_rows; x++)
      {
        mat weights_temp = weights_tmp(x).row(d);
        vec predictions_temp = sum(weights_temp % experts_mat, 1);

        // Sort predictions if quantile_crossing is prohibited
        if (!allow_quantile_crossing)
        {
          predictions_temp = arma::sort(predictions_temp, "ascend", 0);
        }
        predictions_tmp(x).tube(t, d) = predictions_temp;
      }

      // Forecasters prediction
      predictions.tube(t, d) = predictions_tmp(opt_index(t)).tube(t, d);
    }

    for (unsigned int x = 0; x < param_grid.n_rows; x++)
    {

      mat lexp_int(P, K); // Experts loss
      mat lexp_ext(P, K); // Experts loss
      mat lexp(P, K);     // Experts loss
      vec lfor(P);        // Forecasters loss
      cube regret_tmp(D, P, K);
      cube regret(basis_mats_mv(x).n_cols, basis_mats(x).n_cols, K); // Dr x Pr x K

      for (unsigned int d = 0; d < D; d++)
      {
        for (unsigned int p = 0; p < P; p++)
        {
          // Evaluate the ex-post predictive performance
          past_performance(t)(d, p, x) = loss(y(t, d),
                                              predictions_tmp(x)(t - lead_time, d, p),
                                              9999,           // where evaluate loss_gradient
                                              loss_function,  // method
                                              tau(p),         // tau
                                              loss_parameter, // alpha
                                              false);
          if (param_grid(x, 13) != 1)
          {
            for (unsigned int k = 0; k < K; k++)
            {
              lexp_int(p, k) = loss(y(t, d),
                                    experts(t)(d, p, k),
                                    predictions_tmp(x)(t - lead_time, d, p), // where evaluate loss_gradient
                                    loss_function,                           // method
                                    tau(p),                                  // tau
                                    loss_parameter,                          // alpha
                                    loss_gradient);
            }

            if (param_grid(x, 13) == 0)
            {
              lexp.row(p) = lexp_int.row(p);
            }
            else
            {
              lexp_ext.row(p) = arma::vectorise(loss_array.tube(t, p)).t();
              lexp.row(p) = (1 - param_grid(x, 13)) * lexp_int.row(p) + param_grid(x, 13) * lexp_ext.row(p);
            }
          }
          else
          {
            lexp_ext.row(p) = arma::vectorise(loss_array.tube(t, p)).t();
            lexp.row(p) = lexp_ext.row(p);
          }
          lfor(p) = loss(y(t, d),
                         predictions_tmp(x)(t - lead_time, d, p),
                         predictions_tmp(x)(t - lead_time, d, p), // where to evaluate loss_gradient
                         loss_function,                           // method
                         tau(p),                                  // tau
                         loss_parameter,                          // alpha
                         loss_gradient);
        }

        mat regret_int(P, K);
        mat regret_ext(P, K);

        if (param_grid(x, 14) != 1)
        {
          regret_int = (lfor - lexp_int.each_col()).t();
          regret_int *= double(basis_mats(x).n_cols) / double(P);

          if (param_grid(x, 14) == 0)
          {
            regret_tmp.row(d) = regret_int.t();
          }
          else
          {
            regret_ext = regret_array.row(t);
            regret_ext = regret_ext.t();
            regret_ext *= double(basis_mats(x).n_cols) / double(P);
            regret_tmp.row(d) = ((1 - param_grid(x, 14)) * regret_int + param_grid(x, 14) * regret_ext).t();
          }
        }
        else
        {
          regret_ext = regret_array.row(t);
          regret_ext = regret_ext.t();
          regret_ext *= double(basis_mats(x).n_cols) / double(P);
          regret_tmp.row(d) = regret_ext.t();
        }
      }

      for (unsigned int k = 0; k < K; k++)
      {
        mat tmp = (regret_tmp.slice(k) * basis_mats(x)).t() * basis_mats_mv(x);
        regret.slice(k) = tmp.t();
      }

      for (unsigned int dr = 0; dr < regret.n_rows; dr++)
      { // This is subject to change if D will be reduces using another basis

        for (unsigned int pr = 0; pr < regret.n_cols; pr++)
        {

          vec r = regret.tube(dr, pr);

          if (method == "ewa")
          {
            // Update the cumulative regret used by eta
            R(x).tube(dr, pr) = vectorise(R(x).tube(dr, pr) * (1 - param_grid(x, 3))) + r;
            eta(x).tube(dr, pr).fill(param_grid(x, 12));
            beta(x).tube(dr, pr) = vectorise(beta0field(x).tube(dr, pr)).t() * K % softmax_r(param_grid(x, 12) * vectorise(R(x).tube(dr, pr)).t());
          }
          else if (method == "ml_poly")
          {
            // Update the cumulative regret used by ML_Poly
            R(x)
                .tube(dr, pr) = vectorise(R(x).tube(dr, pr) * (1 - param_grid(x, 3))) + r;

            // Update the learning rate
            eta(x).tube(dr, pr) = 1 / (1 / vectorise(eta(x).tube(dr, pr)).t() + square(r.t()));

            // param_grid(x, 12) = gamma
            beta(x).tube(dr, pr) = vectorise(beta0field(x).tube(dr, pr)).t() * K * param_grid(x, 12) % vectorise(eta(x).tube(dr, pr)).t() % pmax_arma(vectorise(R(x).tube(dr, pr)).t(), exp(-700));
            beta(x).tube(dr, pr) /= accu(beta(x).tube(dr, pr));
          }
          else if (method == "boa" || method == "bewa")
          {
            V(x).tube(dr, pr) = vectorise(V(x).tube(dr, pr)).t() * (1 - param_grid(x, 3)) + square(r.t());

            E(x).tube(dr, pr) = max(vectorise(E(x).tube(dr, pr)).t() * (1 - param_grid(x, 3)), abs(r.t()));

            eta(x).tube(dr, pr) =
                pmin_arma(
                    min(1 / (2 * vectorise(E(x).tube(dr, pr))),
                        sqrt(-log(vectorise(beta0field(x).tube(dr, pr))) / vectorise(V(x).tube(dr, pr)))),
                    exp(350));

            vec r_reg = r - vectorise(eta(x).tube(dr, pr)) % square(r);

            R_reg(x).tube(dr, pr) *= (1 - param_grid(x, 3)); // forget
            R_reg(x).tube(dr, pr) +=
                0.5 * (r_reg + conv_to<colvec>::from(vectorise(eta(x).tube(dr, pr)) % r > 0.5) % (2 * vectorise(E(x).tube(dr, pr))));

            if (method == "boa")
            {
              // Wintenberger
              beta(x).tube(dr, pr) = vectorise(beta0field(x).tube(dr, pr)).t() * K % softmax_r(log(param_grid(x, 12) * vectorise(eta(x).tube(dr, pr)).t()) + param_grid(x, 12) * vectorise(eta(x).tube(dr, pr)).t() % vectorise(R_reg(x).tube(dr, pr)).t());
            }
            else
            {
              // Gaillard
              beta(x).tube(dr, pr) = vectorise(beta0field(x).tube(dr, pr)).t() * K % softmax_r(param_grid(x, 12) * vectorise(eta(x).tube(dr, pr)).t() % vectorise(R_reg(x).tube(dr, pr)).t());
            }
          }
          else
          {
            Rcpp::stop("Choose 'boa', 'bewa', 'ml_poly' or 'ewa' as method.");
          }

          //   // Apply thresholds
          if (param_grid(x, 4) > 0)
          {
            int best_k = beta(x).tube(dr, pr).index_max();

            for (double &e : beta(x).tube(dr, pr))
            {
              threshold_soft(e, param_grid(x, 4));
            }
            if (accu(beta(x).tube(dr, pr)) == 0)
            {
              beta(x)(dr, pr, best_k) = 1;
            }
          }

          if (param_grid(x, 5) > 0)
          {
            int best_k = beta(x).tube(dr, pr).index_max();
            for (double &e : beta(x).tube(dr, pr))
            {
              threshold_hard(e, param_grid(x, 5));
            }
            if (accu(beta(x).tube(dr, pr)) == 0)
            {
              beta(x)(dr, pr, best_k) = 1;
            }
          }

          // Add fixed_share
          beta(x).tube(dr, pr) =
              (1 - param_grid(x, 6)) * vectorise(beta(x).tube(dr, pr)) +
              (param_grid(x, 6) / K);
        } // pr
      }   // dr

      cube tmp(D, basis_mats(x).n_cols, K);

      for (unsigned int k = 0; k < K; k++)
      {
        tmp.slice(k) = basis_mats_mv(x) * beta(x).slice(k);

        // // Smoothing
        if (param_grid(x, 7) != -datum::inf)
        {
          // Note that hat was already mutliplied with basis so we can use it directly here
          weights_tmp(x).slice(k) = (hat_mats(x) * tmp.slice(k).t()).t();
        }
        else
        {
          weights_tmp(x).slice(k) = (basis_mats(x) * tmp.slice(k).t()).t();
        }
      }
      // Enshure that constraints hold
      for (unsigned int p = 0; p < P; p++)
      {

        for (unsigned int d = 0; d < D; d++)
        {
          // Positivity
          weights_tmp(x)(span(d), span(p), span::all) =
              pmax_arma(weights_tmp(x)(span(d), span(p), span::all), exp(-700));

          // // Affinity
          weights_tmp(x)(span(d), span(p), span::all) /=
              accu(weights_tmp(x)(span(d), span(p), span::all));
        }
      }
      tmp_performance(x) = accu(past_performance(t).slice(x));
      prog.increment(); // Update progress
      R_CheckUserInterrupt();
    }

    // Sum past_performance in each slice
    cum_performance = (1 - forget_past_performance) * cum_performance + tmp_performance;

    opt_index(t + 1) = cum_performance.index_min();
    chosen_params.row(t) = param_grid.row(opt_index(t + 1));
  }
  // Save Final Weights and Prediction
  weights(T) = weights_tmp(opt_index(T));

  // Predict residual expert forecasts if any are available
  for (unsigned int t = T; t < T + T_E_Y; t++)
  {
    for (unsigned int d = 0; d < D; d++)
    {
      mat experts_mat = experts(t).row(d);
      mat weights_temp = weights(T).row(d);
      vec predictions_temp = sum(weights_temp % experts_mat, 1);

      // Sort predictions if quantile_crossing is prohibited
      if (!allow_quantile_crossing)
      {
        predictions_temp = arma::sort(predictions_temp, "ascend", 0);
      }
      predictions.tube(t, d) = predictions_temp;
    }
  }

  // Save losses suffered by forecaster and experts
  for (unsigned int t = 0; t < T; t++)
  {
    loss_exp(t).set_size(D, P, K);

    for (unsigned int d = 0; d < D; d++)
    {
      for (unsigned int p = 0; p < P; p++)
      {
        for (unsigned int k = 0; k < K; k++)
        {
          loss_exp(t)(d, p, k) =
              loss(y(t, d),
                   experts(t)(d, p, k),
                   9999,           // where to evaluate the loss_gradient
                   loss_function,  // method
                   tau(p),         // tau
                   loss_parameter, // alpha
                   false);         // loss_gradient
        }
        loss_for(t, d, p) = loss(y(t, d),
                                 predictions(t, d, p),
                                 9999,           // where to evaluate the loss_gradient
                                 loss_function,  // method
                                 tau(p),         // tau
                                 loss_parameter, // alpha
                                 false);         // loss_gradient;
      }
    }
  }
  return;
}

// [[Rcpp::export]]
Rcpp::List online_rcpp_mv(
    mat &y,
    arma::field<cube> &experts,
    vec tau, // We don't pass by reference here since tau may be modified
    const unsigned int &lead_time,
    const std::string loss_function,
    const double &loss_parameter,
    const bool &loss_gradient,
    const std::string method,
    const mat &param_grid,
    const double &forget_past_performance,
    bool allow_quantile_crossing,
    const mat w0,
    const mat R0,
    const cube &loss_array,
    const cube &regret_array,
    const bool trace)
{

  // Indexing Convention -> (T, P, K, X)
  // T number of observations
  // D dimensionality of Y. Its 1 in univariate settings
  // P lengths of probability Grid
  // K number of experts
  // X number of parameter combinations to consider

  // Object Dimensions
  const unsigned int T = y.n_rows;
  const unsigned int D = experts(0).n_rows;
  const unsigned int P = experts(0).n_cols;
  const unsigned int K = experts(0).n_slices;
  const unsigned int T_E_Y = experts.n_rows - y.n_rows;

  // Expand y matrix if necessary
  if (y.n_cols == 1)
  {
    y = repmat(y, 1, P);
  }

  // Expand tau if necessary
  if (tau.n_elem == 1)
  {
    tau.resize(P);
    tau.fill(tau(0));
  }

  const unsigned int X = param_grid.n_rows;
  mat chosen_params(T, param_grid.n_cols);
  vec opt_index(T + 1, fill::zeros);
  arma::field<cube> past_performance(T);
  vec tmp_performance(X, fill::zeros);
  vec cum_performance(X, fill::zeros);
  Progress prog(T * X + X, trace);

  // Init objects holding weights
  arma::field<cube> weights_tmp(X);

  // predictions_tmp(T,D P);
  arma::field<cube> predictions_tmp(X);

  // // Output Objects
  cube predictions(T + T_E_Y, D, P, fill::zeros);
  // weights(D, P, K, fill::zeros);
  arma::field<cube> weights(T + 1);

  // // Learning parameters
  arma::field<mat> hat_mats(X);
  arma::field<sp_mat> basis_mats(X);
  arma::field<sp_mat> basis_mats_mv(X);
  arma::field<cube> V(X);
  arma::field<cube> E(X);
  arma::field<cube> k(X);
  arma::field<cube> eta(X);
  arma::field<cube> R(X);
  arma::field<cube> R_reg(X);
  arma::field<cube> beta(X);
  arma::field<cube> beta0field(X);

  vec spline_basis_prob = regspace(1, P) / (P + 1);
  vec spline_basis_mv = regspace(1, D) / (D + 1);

  // Init learning parameters and basis_matrices
  for (unsigned int x = 0; x < X; x++)
  {

    // In second step: skip if recalc is not necessary:
    if (x > 0 &&
        param_grid(x, 2) == param_grid(x - 1, 2) &&
        param_grid(x, 0) == param_grid(x - 1, 0) &&
        param_grid(x, 1) == param_grid(x - 1, 1))
    {
      basis_mats(x) = basis_mats(x - 1);
    }
    else
    {

      basis_mats(x) = make_basis_matrix(spline_basis_prob,
                                        param_grid(x, 0), // kstep
                                        param_grid(x, 2), // degree
                                        param_grid(x, 1), // uneven grid
                                        P % 2 == 0);      // even
    }

    // In second step: skip if recalc is not necessary:
    if (x > 0 &&
        param_grid(x, 15) == param_grid(x - 1, 15) &&
        param_grid(x, 16) == param_grid(x - 1, 16) &&
        param_grid(x, 17) == param_grid(x - 1, 17))
    {
      basis_mats_mv(x) = basis_mats_mv(x - 1);
    }
    else
    {

      basis_mats_mv(x) = make_basis_matrix(spline_basis_mv,
                                           param_grid(x, 15), // kstep
                                           param_grid(x, 16), // degree
                                           param_grid(x, 17), // uneven grid
                                           D % 2 == 0);       // even
    }

    int Pr = basis_mats(x).n_cols;
    int Dr = basis_mats_mv(x).n_cols;

    // Learning parameters
    V(x).zeros(Dr, Pr, K);
    E(x).zeros(Dr, Pr, K);
    k(x).zeros(Dr, Pr, K);

    arma::cube eta_(Dr, Pr, K, fill::zeros);
    eta(x) = eta_;
    if (method == "ml_poly")
    {
      eta_.fill(exp(350));
      eta(x) = eta_;
    }

    R(x).set_size(Dr, Pr, K);
    R_reg(x).set_size(Dr, Pr, K);
    beta(x).set_size(Dr, Pr, K);
    weights_tmp(x).set_size(D, P, K);
    predictions_tmp(x).set_size(T, D, P);

    for (unsigned int dr = 0; dr < Dr; dr++)
    {
      R(x).row(dr) = basis_mats(x).t() * R0;
      R_reg(x).row(dr) = basis_mats(x).t() * R0;
      beta(x).row(dr) = (w0.t() * pinv(mat(basis_mats(x))).t()).t();
    }

    for (unsigned int d = 0; d < D; d++)
    {
      weights_tmp(x).row(d) = w0;
    }
    beta0field(x) = beta(x);

    R_CheckUserInterrupt();
  }

  // Only if smoothing is possible (tau.size > 1)
  if (P > 1)
  {
    // Init hat matrix field
    for (unsigned int x = 0; x < X; x++)
    {
      // In second step: skip if recalc is not necessary:
      if (x > 0 &&
          param_grid(x, 7) == param_grid(x - 1, 7) &&
          param_grid(x, 8) == param_grid(x - 1, 8) &&
          param_grid(x, 10) == param_grid(x - 1, 10) &&
          param_grid(x, 11) == param_grid(x - 1, 11) &&
          param_grid(x, 9) == param_grid(x - 1, 9))
      {
        hat_mats(x) = hat_mats(x - 1);
      }
      else
      {
        if (param_grid(x, 7) != -datum::inf)
          hat_mats(x) = make_hat_matrix(spline_basis_prob,
                                        param_grid(x, 8),  // kstep
                                        param_grid(x, 7),  // lambda
                                        param_grid(x, 11), // differences
                                        param_grid(x, 10), // degree
                                        param_grid(x, 9),  // uneven grid
                                        P % 2 == 0);       // even
      }
      if (param_grid(x, 7) != -datum::inf)
        hat_mats(x) *= basis_mats(x);
      R_CheckUserInterrupt();
      prog.increment(); // Update progress
    }
  }

  // Predictions at t < lead_time using initial weights
  for (unsigned int t = 0; t < lead_time; t++)
  {

    weights(t).set_size(D, P, K);
    past_performance(t).set_size(D, P, X);

    for (unsigned int d = 0; d < D; d++)
    {
      // Save final weights weights_tmp
      weights(t).row(d) = weights_tmp(opt_index(t)).row(d);

      // Store expert predictions temporarily
      mat experts_mat = experts(t).row(d);

      for (unsigned int x = 0; x < X; x++)
      {

        mat weights_temp = weights_tmp(x).row(d);

        // Forecasters prediction
        vec predictions_temp = sum(weights_temp % experts_mat, 1);

        predictions_tmp(x).tube(t, d) = predictions_temp;
      }

      // // Final prediction
      predictions.tube(t, d) = predictions_tmp(opt_index(t)).tube(t, d);
    }

    past_performance(t).fill(datum::nan);
  }

  cube loss_for(T, D, P, fill::zeros);
  field<cube> loss_exp(T);

  const int start = lead_time;

  online_learning_core_mv(
      T,
      D,
      P,
      K,
      T_E_Y,
      start,
      lead_time,
      y,
      experts,
      tau,
      loss_function,
      method,
      loss_parameter,
      loss_gradient,
      weights,
      weights_tmp,
      predictions_tmp,
      predictions,
      past_performance,
      opt_index,
      param_grid,
      chosen_params,
      R,
      R_reg,
      V,
      E,
      eta,
      hat_mats,
      basis_mats,
      basis_mats_mv,
      beta0field,
      beta,
      cum_performance,
      forget_past_performance,
      tmp_performance,
      allow_quantile_crossing,
      loss_for,
      loss_exp,
      loss_array,
      regret_array,
      prog);

  // // 1-Indexing for R-Output
  opt_index = opt_index + 1;

  // Set unused values to NA
  if (lead_time > 0)
  {
    chosen_params.rows(0, lead_time - 1).fill(datum::nan);
  }

  Rcpp::CharacterVector param_names =
      Rcpp::CharacterVector::create("basis_knot_distance",
                                    "basis_knot_distance_power",
                                    "basis_deg",
                                    "forget_regret",
                                    "threshold_soft",
                                    "threshold_hard",
                                    "fixed_share",
                                    "p_smooth_lambda",
                                    "p_smooth_knot_distance",
                                    "p_smooth_knot_distance_power",
                                    "p_smooth_deg",
                                    "smooth_diff",
                                    "gamma",
                                    "loss_share",
                                    "regret_share",
                                    "mv_basis_knot_distance ",
                                    "mv_basis_knot_distance_power",
                                    "mv_basis_deg");

  Rcpp::NumericMatrix parametergrid_out = Rcpp::wrap(param_grid);
  Rcpp::colnames(parametergrid_out) = param_names;

  Rcpp::NumericMatrix chosen_parameters = Rcpp::wrap(chosen_params);
  Rcpp::colnames(chosen_parameters) = param_names;

  Rcpp::List model_data = Rcpp::List::create(
      Rcpp::Named("y") = y,
      Rcpp::Named("experts") = experts,
      Rcpp::Named("tau") = tau);

  Rcpp::List model_parameters = Rcpp::List::create(
      Rcpp::Named("lead_time") = lead_time,
      Rcpp::Named("loss_function") = loss_function,
      Rcpp::Named("loss_parameter") = loss_parameter,
      Rcpp::Named("loss_gradient") = loss_gradient,
      Rcpp::Named("method") = method,
      Rcpp::Named("forget_past_performance") = forget_past_performance,
      Rcpp::Named("allow_quantile_crossing") = allow_quantile_crossing);

  Rcpp::List model_objects = Rcpp::List::create(
      Rcpp::Named("weights_tmp") = weights_tmp,
      Rcpp::Named("predictions_tmp") = predictions_tmp,
      Rcpp::Named("cum_performance") = cum_performance,
      Rcpp::Named("hat_matrices") = hat_mats,
      Rcpp::Named("V") = V,
      Rcpp::Named("E") = E,
      Rcpp::Named("k") = k,
      Rcpp::Named("eta") = eta,
      Rcpp::Named("R") = R,
      Rcpp::Named("R_reg") = R_reg,
      Rcpp::Named("beta") = beta,
      Rcpp::Named("beta0field") = beta0field);

  Rcpp::List model_spec = Rcpp::List::create(
      Rcpp::Named("data") = model_data,
      Rcpp::Named("parameters") = model_parameters,
      Rcpp::Named("objects") = model_objects);

  Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("predictions") = predictions,
      Rcpp::Named("weights") = weights,
      Rcpp::Named("forecaster_loss") = loss_for,
      Rcpp::Named("experts_loss") = loss_exp,
      Rcpp::Named("past_performance") = past_performance,
      Rcpp::Named("chosen_parameters") = chosen_parameters,
      Rcpp::Named("opt_index") = opt_index,
      Rcpp::Named("parametergrid") = parametergrid_out,
      Rcpp::Named("basis_matrices") = basis_mats,
      Rcpp::Named("specification") = model_spec);

  // Rcpp::List out;
  out.attr("class") = "online";

  return out;
}
