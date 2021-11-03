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
    const unsigned int &P,
    const unsigned int &K,
    const unsigned int &T_E_Y,
    const unsigned int &start,
    const unsigned int &lead_time,
    const mat &y,
    const cube &experts,
    const vec &tau,
    const std::string loss_function,
    const std::string method,
    const double &loss_parameter,
    const bool &loss_gradient,
    cube &weights,
    cube &weights_tmp,
    cube &predictions_tmp,
    mat &predictions,
    cube &past_performance,
    vec &opt_index,
    const mat &param_grid,
    mat &chosen_params,
    field<mat> &R,
    field<mat> &R_reg,
    field<mat> &V,
    field<mat> &E,
    field<mat> &eta,
    field<mat> &hat_mats,
    field<sp_mat> &basis_mats,
    field<mat> &beta0field,
    field<mat> &beta,
    vec &cum_performance,
    const double &forget_past_performance,
    vec &tmp_performance,
    const bool &allow_quantile_crossing,
    mat &loss_for,
    cube &loss_exp,
    const cube &loss_array,
    const cube &regret_array,
    Progress &prog)
{
  for (unsigned int t = start; t < T; t++)
  {

    // Save final weights weights_tmp
    weights.row(t) = weights_tmp.slice(opt_index(t));

    // Store expert predictions temporarily
    mat experts_mat = experts.row(t);

    // Forecasters prediction
    mat predictions_temp = sum(weights_tmp.each_slice() % experts_mat, 1);

    // Sort predictions if quantile_crossing is prohibited
    if (!allow_quantile_crossing)
    {
      predictions_temp = arma::sort(predictions_temp, "ascend", 0);
    }

    predictions_tmp.row(t) = predictions_temp;

    // Final prediction
    predictions.row(t) =
        vectorise(predictions_tmp(span(t), span::all, span(opt_index(t)))).t();

    for (unsigned int x = 0; x < param_grid.n_rows; x++)
    {

      mat lexp_int(P, K); // Experts loss
      mat lexp_ext(P, K); // Experts loss
      mat lexp(P, K);     // Experts loss
      vec lfor(P);        // Forecasters loss

      for (unsigned int p = 0; p < P; p++)
      {

        // Evaluate the ex-post predictive performance
        past_performance(t, p, x) = loss(y(t, p),
                                         predictions_tmp(t - lead_time, p, x),
                                         9999,           // where evaluate loss_gradient
                                         loss_function,  // method
                                         tau(p),         // tau
                                         loss_parameter, // alpha
                                         false);

        if (param_grid(x, 13) != 1)
        {
          for (unsigned int k = 0; k < K; k++)
          {
            lexp_int(p, k) = loss(y(t, p),
                                  experts(t, p, k),
                                  predictions_tmp(t - lead_time, p, x), // where evaluate loss_gradient
                                  loss_function,                        // method
                                  tau(p),                               // tau
                                  loss_parameter,                       // alpha
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

        lfor(p) = loss(y(t, p),
                       predictions_tmp(t - lead_time, p, x),
                       predictions_tmp(t - lead_time, p, x), // where to evaluate loss_gradient
                       loss_function,                        // method
                       tau(p),                               // tau
                       loss_parameter,                       // alpha
                       loss_gradient);
      }

      mat regret_int;
      mat regret_ext;
      mat regret;

      if (param_grid(x, 14) != 1)
      {
        regret_int = (lfor - lexp_int.each_col()).t();
        regret_int *= double(basis_mats(x).n_cols) / double(P);

        if (param_grid(x, 14) == 0)
        {
          regret = regret_int;
        }
        else
        {
          regret_ext = regret_array.row(t);
          regret_ext = regret_ext.t();
          regret_ext *= double(basis_mats(x).n_cols) / double(P);
          regret = (1 - param_grid(x, 14)) * regret_int + param_grid(x, 14) * regret_ext;
        }
      }
      else
      {
        regret_ext = regret_array.row(t);
        regret_ext = regret_ext.t();
        regret_ext *= double(basis_mats(x).n_cols) / double(P);
        regret = regret_ext;
      }

      regret *= basis_mats(x);

      for (unsigned int l = 0; l < regret.n_cols; l++)
      {

        vec r = regret.col(l);

        if (method == "ewa")
        {
          // Update the cumulative regret used by eta
          R(x).row(l) = R(x).row(l) * (1 - param_grid(x, 3)) + r.t();
          eta(x).row(l).fill(param_grid(x, 12));
          beta(x).row(l) = beta0field(x).row(l) * K % softmax_r(param_grid(x, 12) * R(x).row(l));
        }
        else if (method == "ml_poly")
        {
          // Update the cumulative regret used by ML_Poly
          R(x).row(l) = R(x).row(l) * (1 - param_grid(x, 3)) + r.t();

          // Update the learning rate
          eta(x).row(l) = 1 / (1 / eta(x).row(l) + square(r.t()));

          // param_grid(x, 12) = gamma
          beta(x).row(l) =
              beta0field(x).row(l) * K * param_grid(x, 12) % eta(x).row(l) % pmax_arma(R(x).row(l), exp(-700));
          beta(x).row(l) /= accu(beta(x).row(l));
        }
        else if (method == "boa" || method == "bewa")
        {

          V(x).row(l) = V(x).row(l) * (1 - param_grid(x, 3)) + square(r.t());

          E(x).row(l) = max(E(x).row(l) * (1 - param_grid(x, 3)), abs(r.t()));

          eta(x).row(l) =
              pmin_arma(
                  min(1 / (2 * E(x).row(l)),
                      sqrt(-log(beta0field(x).row(l)) / V(x).row(l))),
                  exp(350));

          vec r_reg = r - eta(x).row(l).t() % square(r);

          R_reg(x).row(l) *= (1 - param_grid(x, 3));
          R_reg(x).row(l) +=
              0.5 * (r_reg.t() + conv_to<colvec>::from(eta(x).row(l) % r.t() > 0.5).t() % (2 * E(x).row(l)));

          if (method == "boa")
          {
            // Wintenberger
            beta(x).row(l) = beta0field(x).row(l) * K % softmax_r(log(param_grid(x, 12) * eta(x).row(l)) + param_grid(x, 12) * eta(x).row(l) % R_reg(x).row(l));
          }
          else
          {
            // Gaillard
            beta(x).row(l) = beta0field(x).row(l) * K % softmax_r(param_grid(x, 12) * eta(x).row(l) % R_reg(x).row(l));
          }
        }
        else
        {
          Rcpp::stop("Choose 'boa', 'bewa', 'ml_poly' or 'ewa' as method.");
        }

        // Apply thresholds
        if (param_grid(x, 4) > 0)
        {
          int best_k = beta(x).row(l).index_max();
          for (double &e : beta(x).row(l))
          {
            threshold_soft(e, param_grid(x, 4));
          }
          if (accu(beta(x).row(l)) == 0)
          {
            beta(x)(l, best_k) = 1;
          }
        }

        if (param_grid(x, 5) > 0)
        {
          int best_k = beta(x).row(l).index_max();
          for (double &e : beta(x).row(l))
          {
            threshold_hard(e, param_grid(x, 5));
          }
          if (accu(beta(x).row(l)) == 0)
          {
            beta(x)(l, best_k) = 1;
          }
        }

        // Add fixed_share
        beta(x).row(l) =
            (1 - param_grid(x, 6)) * beta(x).row(l) +
            (param_grid(x, 6) / K);
      }

      // Smoothing
      if (param_grid(x, 7) != -datum::inf)
      {
        // Note that hat was already mutliplied with basis so we can use it directly here
        weights_tmp.slice(x) = hat_mats(x) * beta(x);
      }
      else
      {
        weights_tmp.slice(x) = basis_mats(x) * beta(x);
      }

      // Enshure that constraints hold
      for (unsigned int p = 0; p < P; p++)
      {

        // Positivity
        weights_tmp(span(p), span::all, span(x)) =
            pmax_arma(weights_tmp(span(p), span::all, span(x)), exp(-700));

        // Affinity
        weights_tmp(span(p), span::all, span(x)) /=
            accu(weights_tmp(span(p), span::all, span(x)));
      }

      tmp_performance(x) = accu(past_performance(span(t), span::all, span(x)));
      prog.increment(); // Update progress
      R_CheckUserInterrupt();
    }

    // Sum past_performance in each slice
    cum_performance = (1 - forget_past_performance) * cum_performance + tmp_performance;
    opt_index(t + 1) = cum_performance.index_min();
    chosen_params.row(t) = param_grid.row(opt_index(t + 1));
  }

  // Save Weights and Prediction
  weights.row(T) = weights_tmp.slice(opt_index(T));

  // Predict residual expert forecasts if any are available
  for (unsigned int t = T; t < T + T_E_Y; t++)
  {
    mat experts_mat = experts.row(t);
    mat predictions_temp = sum(weights_tmp.slice(opt_index(T)) % experts_mat, 1);
    // Sort predictions if quantile_crossing is prohibited
    if (!allow_quantile_crossing)
    {
      predictions_temp = arma::sort(predictions_temp, "ascend", 0);
    }
    predictions.row(t) = vectorise(predictions_temp).t();
  }

  // Save losses suffered by forecaster and experts
  for (unsigned int t = 0; t < T; t++)
  {
    for (unsigned int p = 0; p < P; p++)
    {
      for (unsigned int k = 0; k < K; k++)
      {
        loss_exp(t, p, k) =
            loss(y(t, p),
                 experts(t, p, k),
                 9999,           // where to evaluate the loss_gradient
                 loss_function,  // method
                 tau(p),         // tau
                 loss_parameter, // alpha
                 false);         // loss_gradient
      }
      loss_for(t, p) = loss(y(t, p),
                            predictions(t, p),
                            9999,           // where to evaluate the loss_gradient
                            loss_function,  // method
                            tau(p),         // tau
                            loss_parameter, // alpha
                            false);         // loss_gradient;
    }
  }

  return;
}

// [[Rcpp::export]]
Rcpp::List online_rcpp_mv(
    mat &y,
    field<cube> &experts,
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
  // if (y.n_cols == 1)
  // {
  //   y = repmat(y, 1, P);
  // }

  // // Expand tau if necessary
  // if (tau.n_elem == 1)
  // {
  //   tau.resize(P);
  //   tau.fill(tau(0));
  // }

  // const unsigned int X = param_grid.n_rows;
  // mat chosen_params(T, param_grid.n_cols);
  // vec opt_index(T + 1, fill::zeros);
  // cube past_performance(T, P, X, fill::zeros);
  // vec tmp_performance(X, fill::zeros);
  // vec cum_performance(X, fill::zeros);
  // Progress prog(T * X + X, trace);

  // // Init objects holding weights
  // cube weights_tmp(P, K, X);
  // weights_tmp.each_slice() = w0;

  // cube predictions_tmp(T, P, X);

  // // Output Objects
  // mat predictions(T + T_E_Y, P, fill::zeros);
  // cube weights(T + 1, P, K, fill::zeros);

  // // Learning parameters
  // field<mat> hat_mats(X);
  // field<sp_mat> basis_mats(X);
  // field<mat> V(X);
  // field<mat> E(X);
  // field<mat> k(X);
  // field<mat> eta(X);
  // field<mat> R(X);
  // field<mat> R_reg(X);
  // field<mat> beta(X);
  // field<mat> beta0field(X);

  // vec spline_basis_x = regspace(1, P) / (P + 1);

  // // Init hat matrix field
  // for (unsigned int x = 0; x < X; x++)
  // {

  //   // In second step: skip if recalc is not necessary:
  //   if (x > 0 &&
  //       param_grid(x, 2) == param_grid(x - 1, 2) &&
  //       param_grid(x, 0) == param_grid(x - 1, 0) &&
  //       param_grid(x, 1) == param_grid(x - 1, 1))
  //   {
  //     basis_mats(x) = basis_mats(x - 1);
  //   }
  //   else
  //   {

  //     basis_mats(x) = make_basis_matrix(spline_basis_x,
  //                                       param_grid(x, 0), // kstep
  //                                       param_grid(x, 2), // degree
  //                                       param_grid(x, 1), // uneven grid
  //                                       P % 2 == 0);      // even
  //   }

  //   int L = basis_mats(x).n_cols;
  //   V(x).zeros(L, K);
  //   E(x).zeros(L, K);
  //   k(x).zeros(L, K);

  //   mat eta_(L, K, fill::zeros);
  //   eta(x) = eta_;
  //   if (method == "ml_poly")
  //   {
  //     eta_.fill(exp(350));
  //     eta(x) = eta_;
  //   }

  //   R_reg(x) = basis_mats(x).t() * R0;
  //   R(x) = basis_mats(x).t() * R0;

  //   beta(x) = (w0.t() * pinv(mat(basis_mats(x))).t()).t();

  //   beta0field(x) = beta(x);

  //   R_CheckUserInterrupt();
  // }

  // // Only if smoothing is possible (tau.size > 1)
  // if (P > 1)
  // {
  //   // Init hat matrix field
  //   for (unsigned int x = 0; x < X; x++)
  //   {
  //     // In second step: skip if recalc is not necessary:
  //     if (x > 0 &&
  //         param_grid(x, 7) == param_grid(x - 1, 7) &&
  //         param_grid(x, 8) == param_grid(x - 1, 8) &&
  //         param_grid(x, 10) == param_grid(x - 1, 10) &&
  //         param_grid(x, 11) == param_grid(x - 1, 11) &&
  //         param_grid(x, 9) == param_grid(x - 1, 9))
  //     {
  //       hat_mats(x) = hat_mats(x - 1);
  //     }
  //     else
  //     {
  //       if (param_grid(x, 7) != -datum::inf)
  //         hat_mats(x) = make_hat_matrix(spline_basis_x,
  //                                       param_grid(x, 8),  // kstep
  //                                       param_grid(x, 7),  // lambda
  //                                       param_grid(x, 11), // differences
  //                                       param_grid(x, 10), // degree
  //                                       param_grid(x, 9),  // uneven grid
  //                                       P % 2 == 0);       // even
  //     }
  //     if (param_grid(x, 7) != -datum::inf)
  //       hat_mats(x) *= basis_mats(x);
  //     R_CheckUserInterrupt();
  //     prog.increment(); // Update progress
  //   }
  // }

  // // Predictions at t < lead_time using initial weights
  // for (unsigned int t = 0; t < lead_time; t++)
  // {
  //   // Save final weights weights_tmp
  //   weights.row(t) = weights_tmp.slice(opt_index(t));

  //   // Store expert predictions temporarily
  //   mat experts_mat = experts.row(t);

  //   // Forecasters prediction
  //   mat predictions_temp = sum(weights_tmp.each_slice() % experts_mat, 1);
  //   predictions_tmp.row(t) = predictions_temp;

  //   // Final prediction
  //   predictions.row(t) =
  //       vectorise(predictions_tmp(span(t), span::all, span(opt_index(t)))).t();

  //   past_performance.row(t).fill(datum::nan);
  // }

  // mat loss_for(T, P, fill::zeros);
  // cube loss_exp(T, P, K, fill::zeros);

  // const int start = lead_time;

  // online_learning_core(
  //     T,
  //     P,
  //     K,
  //     T_E_Y,
  //     start,
  //     lead_time,
  //     y,
  //     experts,
  //     tau,
  //     loss_function,
  //     method,
  //     loss_parameter,
  //     loss_gradient,
  //     weights,
  //     weights_tmp,
  //     predictions_tmp,
  //     predictions,
  //     past_performance,
  //     opt_index,
  //     param_grid,
  //     chosen_params,
  //     R,
  //     R_reg,
  //     V,
  //     E,
  //     eta,
  //     hat_mats,
  //     basis_mats,
  //     beta0field,
  //     beta,
  //     cum_performance,
  //     forget_past_performance,
  //     tmp_performance,
  //     allow_quantile_crossing,
  //     loss_for,
  //     loss_exp,
  //     loss_array,
  //     regret_array,
  //     prog);

  // // 1-Indexing for R-Output
  // opt_index = opt_index + 1;

  // // Set unused values to NA
  // if (lead_time > 0)
  // {
  //   chosen_params.rows(0, lead_time - 1).fill(datum::nan);
  // }

  // Rcpp::CharacterVector param_names =
  //     Rcpp::CharacterVector::create("basis_knot_distance",
  //                                   "basis_knot_distance_power",
  //                                   "basis_deg",
  //                                   "forget_regret",
  //                                   "threshold_soft",
  //                                   "threshold_hard",
  //                                   "fixed_share",
  //                                   "p_smooth_lambda",
  //                                   "p_smooth_knot_distance",
  //                                   "p_smooth_knot_distance_power",
  //                                   "p_smooth_deg",
  //                                   "smooth_diff",
  //                                   "gamma",
  //                                   "loss_share",
  //                                   "regret_share");

  // Rcpp::NumericMatrix parametergrid_out = Rcpp::wrap(param_grid);
  // Rcpp::colnames(parametergrid_out) = param_names;

  // Rcpp::NumericMatrix chosen_parameters = Rcpp::wrap(chosen_params);
  // Rcpp::colnames(chosen_parameters) = param_names;

  // Rcpp::List model_data = Rcpp::List::create(
  //     Rcpp::Named("y") = y,
  //     Rcpp::Named("experts") = experts,
  //     Rcpp::Named("tau") = tau);

  // Rcpp::List model_parameters = Rcpp::List::create(
  //     Rcpp::Named("lead_time") = lead_time,
  //     Rcpp::Named("loss_function") = loss_function,
  //     Rcpp::Named("loss_parameter") = loss_parameter,
  //     Rcpp::Named("loss_gradient") = loss_gradient,
  //     Rcpp::Named("method") = method,
  //     Rcpp::Named("forget_past_performance") = forget_past_performance,
  //     Rcpp::Named("allow_quantile_crossing") = allow_quantile_crossing);

  // Rcpp::List model_objects = Rcpp::List::create(
  //     Rcpp::Named("weights_tmp") = weights_tmp,
  //     Rcpp::Named("predictions_tmp") = predictions_tmp,
  //     Rcpp::Named("cum_performance") = cum_performance,
  //     Rcpp::Named("hat_matrices") = hat_mats,
  //     Rcpp::Named("V") = V,
  //     Rcpp::Named("E") = E,
  //     Rcpp::Named("k") = k,
  //     Rcpp::Named("eta") = eta,
  //     Rcpp::Named("R") = R,
  //     Rcpp::Named("R_reg") = R_reg,
  //     Rcpp::Named("beta") = beta,
  //     Rcpp::Named("beta0field") = beta0field);

  // Rcpp::List model_spec = Rcpp::List::create(
  //     Rcpp::Named("data") = model_data,
  //     Rcpp::Named("parameters") = model_parameters,
  //     Rcpp::Named("objects") = model_objects);

  // Rcpp::List out = Rcpp::List::create(
  //     Rcpp::Named("predictions") = predictions,
  //     Rcpp::Named("weights") = weights,
  //     Rcpp::Named("forecaster_loss") = loss_for,
  //     Rcpp::Named("experts_loss") = loss_exp,
  //     Rcpp::Named("past_performance") = past_performance,
  //     Rcpp::Named("chosen_parameters") = chosen_parameters,
  //     Rcpp::Named("opt_index") = opt_index,
  //     Rcpp::Named("parametergrid") = parametergrid_out,
  //     Rcpp::Named("basis_matrices") = basis_mats,
  //     Rcpp::Named("specification") = model_spec);

  // out.attr("class") = "online";
  Rcpp::List out;

  return out;
}
