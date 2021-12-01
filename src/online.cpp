// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppClock)]]
// [[Rcpp::plugins("cpp11")]]

#include <misc.h>
#include <loss.h>
#include <splines2.h>
#include <splines.h>

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <RcppClock.h>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

void online_learning_core(
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
    field<mat> &hat_mats_mv,
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
    const field<cube> &loss_array,
    const field<cube> &regret_array,
    Progress &prog,
    Rcpp::Clock &clock)
{
  for (unsigned int t = start; t < T; t++)
  {

    weights(t).set_size(D, P, K);
    past_performance(t).set_size(D, P, param_grid.n_rows);

    clock.tick("loss");
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

    clock.tock("loss");
    for (unsigned int x = 0; x < param_grid.n_rows; x++)
    {
      clock.tick("regret");
      mat lexp_int(P, K); // Experts loss
      mat lexp_ext(P, K); // Experts loss
      mat lexp(P, K);     // Experts loss
      vec lfor(P);        // Forecasters loss
      cube regret_tmp(D, P, K);
      cube regret(basis_mats_mv(x).n_rows, basis_mats(x).n_cols, K); // Dr x Pr x K

      for (unsigned int d = 0; d < D; d++)
      {
#pragma omp parallel for
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
              lexp_ext.row(p) = arma::vectorise(loss_array(t).tube(d, p)).t();
              lexp.row(p) = (1 - param_grid(x, 13)) * lexp_int.row(p) + param_grid(x, 13) * lexp_ext.row(p);
            }
          }
          else
          {
            lexp_ext.row(p) = arma::vectorise(loss_array(t).tube(d, p)).t();
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
            regret_ext = regret_array(t).row(d);
            regret_ext = regret_ext.t();
            regret_ext *= double(basis_mats(x).n_cols) / double(P);
            regret_tmp.row(d) = ((1 - param_grid(x, 14)) * regret_int + param_grid(x, 14) * regret_ext).t();
          }
        }
        else
        {
          regret_ext = regret_array(t).row(d);
          regret_ext = regret_ext.t();
          regret_ext *= double(basis_mats(x).n_cols) / double(P);
          regret_tmp.row(d) = regret_ext.t();
        }
      }

#pragma omp parallel for
      for (unsigned int k = 0; k < K; k++)
      {
        regret.slice(k) = basis_mats_mv(x) * regret_tmp.slice(k) * basis_mats(x);
      }

      clock.tock("regret");
      clock.tick("learning");
#pragma omp parallel for collapse(2)
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
      clock.tock("learning");

      clock.tick("smoothing");
#pragma omp parallel for
      // Smoothing
      for (unsigned int k = 0; k < K; k++)
      {
        if (param_grid(x, 7) != -datum::inf && param_grid(x, 18) != -datum::inf)
        {
          weights_tmp(x).slice(k) = hat_mats_mv(x) * beta(x).slice(k) * hat_mats(x);
        }
        else if (param_grid(x, 7) != -datum::inf)
        {
          weights_tmp(x).slice(k) = basis_mats_mv(x).t() * beta(x).slice(k) * hat_mats(x);
        }
        else if (param_grid(x, 18) != -datum::inf)
        {
          weights_tmp(x).slice(k) = hat_mats_mv(x) * beta(x).slice(k) * basis_mats(x).t();
        }
        else
        {
          weights_tmp(x).slice(k) = basis_mats_mv(x).t() * beta(x).slice(k) * basis_mats(x).t();
        }
      }
      clock.tock("smoothing");

#pragma omp parallel for
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
      R_CheckUserInterrupt();
    }

    // Sum past_performance in each slice
    cum_performance = (1 - forget_past_performance) * cum_performance + tmp_performance;

    opt_index(t + 1) = cum_performance.index_min();
    chosen_params.row(t) = param_grid.row(opt_index(t + 1));
    prog.increment(); // Update progress

  } // t

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
  clock.tick("loss_for_exp");
#pragma omp parallel for
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
  clock.tock("loss_for_exp");
  return;
}

// [[Rcpp::export]]
Rcpp::List online_rcpp(
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
    const cube w0,
    const cube R0,
    const field<cube> &loss_array,
    const field<cube> &regret_array,
    const bool trace)
{

  Rcpp::Clock clock;
  clock.tick("init");

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
  Progress prog(T + 2 * X, trace);

  // Init objects holding weights
  arma::field<cube> weights_tmp(X);
  arma::field<cube> predictions_tmp(X); // TODO Change to T as first dim

  // Output Objects
  cube predictions(T + T_E_Y, D, P, fill::zeros);
  arma::field<cube> weights(T + 1);

  // Learning parameters
  arma::field<mat> hat_mats(X);
  arma::field<mat> hat_mats_mv(X);
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

  // Init hat and basis_matrices
#pragma omp parallel for
  for (unsigned int x = 0; x < X; x++)
  {

    // In second step: skip if recalc is not necessary:
    if (x == 0)
    {
      basis_mats(x) = make_basis_matrix(spline_basis_prob,
                                        param_grid(x, 0), // kstep
                                        param_grid(x, 2), // degree
                                        param_grid(x, 1), // uneven grid
                                        P % 2 == 0);      // even
    }
    else if (
        param_grid(x, 2) != param_grid(x - 1, 2) ||
        param_grid(x, 0) != param_grid(x - 1, 0) ||
        param_grid(x, 1) != param_grid(x - 1, 1))
    {
      basis_mats(x) = make_basis_matrix(spline_basis_prob,
                                        param_grid(x, 0), // kstep
                                        param_grid(x, 2), // degree
                                        param_grid(x, 1), // uneven grid
                                        P % 2 == 0);      // even
    }

    // In second step: skip if recalc is not necessary:
    if (x == 0)
    {
      basis_mats_mv(x) = make_basis_matrix(spline_basis_mv,
                                           param_grid(x, 15), // kstep
                                           param_grid(x, 16), // degree
                                           param_grid(x, 17), // uneven grid
                                           D % 2 == 0)
                             .t(); // even
    }
    else if (param_grid(x, 15) != param_grid(x - 1, 15) ||
             param_grid(x, 16) != param_grid(x - 1, 16) ||
             param_grid(x, 17) != param_grid(x - 1, 17))
    {
      basis_mats_mv(x) = make_basis_matrix(spline_basis_mv,
                                           param_grid(x, 15), // kstep
                                           param_grid(x, 16), // degree
                                           param_grid(x, 17), // uneven grid
                                           D % 2 == 0)
                             .t(); // even
    }

    if (x == 0)
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
    else if (
        param_grid(x, 7) != param_grid(x - 1, 7) ||
        param_grid(x, 8) != param_grid(x - 1, 8) ||
        param_grid(x, 10) != param_grid(x - 1, 10) ||
        param_grid(x, 11) != param_grid(x - 1, 11) ||
        param_grid(x, 9) != param_grid(x - 1, 9))
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

    if (x == 0 ||
        param_grid(x, 18) != param_grid(x - 1, 18) ||
        param_grid(x, 19) != param_grid(x - 1, 19) ||
        param_grid(x, 20) != param_grid(x - 1, 20) ||
        param_grid(x, 21) != param_grid(x - 1, 21) ||
        param_grid(x, 22) != param_grid(x - 1, 22))
    {
      if (param_grid(x, 18) != -datum::inf)
        hat_mats_mv(x) = make_hat_matrix(spline_basis_mv,
                                         param_grid(x, 19), // kstep
                                         param_grid(x, 18), // lambda
                                         param_grid(x, 22), // differences
                                         param_grid(x, 21), // degree
                                         param_grid(x, 20), // uneven grid
                                         P % 2 == 0);       // even
    }
    prog.increment();
  }

  // Note that this can't be parralelized due to basis_mats(x - 1)
  for (unsigned int x = 0; x < X; x++)
  {
    if (x > 0)
    {
      if (
          param_grid(x, 2) == param_grid(x - 1, 2) &&
          param_grid(x, 0) == param_grid(x - 1, 0) &&
          param_grid(x, 1) == param_grid(x - 1, 1))
      {
        basis_mats(x) = basis_mats(x - 1);
      }

      if (
          param_grid(x, 7) == param_grid(x - 1, 7) &&
          param_grid(x, 8) == param_grid(x - 1, 8) &&
          param_grid(x, 10) == param_grid(x - 1, 10) &&
          param_grid(x, 11) == param_grid(x - 1, 11) &&
          param_grid(x, 9) == param_grid(x - 1, 9))
      {
        hat_mats(x) = hat_mats(x - 1);
      }

      if (
          param_grid(x, 15) == param_grid(x - 1, 15) &&
          param_grid(x, 16) == param_grid(x - 1, 16) &&
          param_grid(x, 17) == param_grid(x - 1, 17))
      {
        basis_mats_mv(x) = basis_mats_mv(x - 1);
      }

      if (
          param_grid(x, 18) == param_grid(x - 1, 18) &&
          param_grid(x, 19) == param_grid(x - 1, 19) &&
          param_grid(x, 20) == param_grid(x - 1, 20) &&
          param_grid(x, 21) == param_grid(x - 1, 21) &&
          param_grid(x, 22) == param_grid(x - 1, 22))
      {
        hat_mats_mv(x) = hat_mats_mv(x - 1);
      }
    }
  }

#pragma omp parallel for
  for (unsigned int x = 0; x < X; x++)
  {

    if (param_grid(x, 7) != -datum::inf)
    {
      hat_mats(x) = basis_mats(x).t() * hat_mats(x);
    }

    if (param_grid(x, 18) != -datum::inf)
    {
      hat_mats_mv(x) = hat_mats_mv(x) * basis_mats_mv(x).t();
    }

    unsigned int Pr = basis_mats(x).n_cols;
    unsigned int Dr = basis_mats_mv(x).n_rows;

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

    for (unsigned int d = 0; d < D; d++)
    {
      weights_tmp(x).row(d) = w0.row(d);
    }

    for (unsigned int k = 0; k < K; k++)
    {
      R(x).slice(k) = basis_mats_mv(x) * R0.slice(k) * basis_mats(x);
      R_reg(x).slice(k) = basis_mats_mv(x) * R0.slice(k) * basis_mats(x);
      beta(x).slice(k) = pinv(mat(basis_mats_mv(x))).t() * w0.slice(k) * pinv(mat(basis_mats(x))).t();
    }

    beta0field(x) = beta(x);
    prog.increment(); // Update progress
  }

  R_CheckUserInterrupt();

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
  clock.tock("init");

  clock.tick("core");

  online_learning_core(
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
      hat_mats_mv,
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
      prog,
      clock);

  clock.tock("core");

  clock.tick("wrangle");

  // // 1-Indexing for R-Output
  opt_index = opt_index + 1;

  // Set unused values to NA
  if (lead_time > 0)
  {
    chosen_params.rows(0, lead_time - 1).fill(datum::nan);
  }

  Rcpp::NumericMatrix parametergrid_out = Rcpp::wrap(param_grid);

  Rcpp::NumericMatrix chosen_parameters = Rcpp::wrap(chosen_params);

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
      Rcpp::Named("hat_matrices_mv") = hat_mats_mv,
      Rcpp::Named("basis_matrices") = basis_mats,
      Rcpp::Named("basis_matrices_mv") = basis_mats_mv,
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
      Rcpp::Named("specification") = model_spec);

  // Rcpp::List out;
  out.attr("class") = "online";

  clock.tock("wrangle");

  clock.stop("times");

  // Rcpp::List out;
  return out;
}

// [[Rcpp::export]]
Rcpp::List predict_online(
    Rcpp::List &object,
    cube &new_experts)
{

  // This creates a reference to the specification, not a copy
  Rcpp::List specification = object["specification"];

  Rcpp::List model_parameters = specification["parameters"];
  Rcpp::List model_data = specification["data"];

  bool allow_quantile_crossing = model_parameters["allow_quantile_crossing"];
  cube experts = model_data["experts"];
  experts.insert_rows(experts.n_rows, new_experts);
  model_data["experts"] = experts;

  mat predictions = object["predictions"];
  cube weights = object["weights"];

  mat predictions_new(new_experts.n_rows, new_experts.n_cols);

  // Sort predictions if quantile_crossing is prohibited
  if (!allow_quantile_crossing)
  {
    predictions_new = arma::sort(predictions_new, "ascend", 1);
  }

  mat predictions_joined = join_vert(predictions, predictions_new);

  object["predictions"] = predictions_joined;

  return object;
}

// [[Rcpp::export]]
Rcpp::List update_online(
    Rcpp::List &object,
    mat &new_y,
    arma::field<cube> &new_experts)
{
  Rcpp::Clock clock;

  clock.tick("init");
  const unsigned int new_T = new_experts.n_rows;
  // This creates a references not copies
  Rcpp::List specification = object["specification"];
  Rcpp::List model_parameters = specification["parameters"];
  Rcpp::List model_data = specification["data"];
  Rcpp::List model_objects = specification["objects"];

  // Data
  arma::field<cube> experts = model_data["experts"];
  const unsigned int old_T = experts.n_rows;

  // Join old and new expert_predictions
  arma::field<cube> experts_new(old_T + new_T);
  if (new_T > 0)
  {
    experts_new.rows(0, old_T - 1) = experts;
    experts_new.rows(old_T, old_T + new_T - 1) = new_experts;
  }
  else
  {
    experts_new = experts;
  }

  mat y = model_data["y"];
  y.insert_rows(y.n_rows, new_y);

  Progress prog(999, false); // TODO

  // Object Dimensions
  const unsigned int T = y.n_rows;
  const unsigned int D = experts(0).n_rows;
  const unsigned int P = experts(0).n_cols;
  const unsigned int K = experts(0).n_slices;
  const int T_E_Y = experts_new.n_rows - y.n_rows;

  if (T_E_Y < 0)
    Rcpp::stop("Number of provided expert predictions has to match or exceed observations.");

  vec tau = model_data["tau"];
  mat param_grid = object["parametergrid"];
  const int X = param_grid.n_rows;
  mat chosen_params = object["chosen_parameters"];
  chosen_params.resize(T, param_grid.n_cols);
  vec opt_index = object["opt_index"];
  // Zero indexing in C++
  opt_index -= 1;
  opt_index.resize(T + 1);

  arma::field<cube> past_performance(T);
  past_performance.rows(0, old_T - 1) =
      Rcpp::as<arma::field<cube>>(object["past_performance"]);

  vec tmp_performance(X, fill::zeros);
  vec cum_performance = model_objects["cum_performance"];

  arma::field<cube> weights_tmp =
      Rcpp::as<arma::field<cube>>(model_objects["weights_tmp"]);

  arma::field<cube> predictions_tmp =
      Rcpp::as<arma::field<cube>>(model_objects["predictions_tmp"]);

  for (unsigned int x = 0; x < X; x++)
  {
    predictions_tmp(x).resize(T + T_E_Y, D, P);
  }

  // TODO Add loss_array and regret_array functionality

  arma::field<cube> loss_array;
  arma::field<cube> regret_array;

  // Output Objects
  cube predictions = object["predictions"];
  predictions.resize(T + T_E_Y, D, P);
  arma::field<cube> weights(T + 1);
  weights.rows(0, old_T) = Rcpp::as<arma::field<cube>>(object["weights"]);

  field<mat> hat_mats = model_objects["hat_matrices"];
  field<sp_mat> basis_mats = model_objects["basis_matrices"];
  field<mat> hat_mats_mv = model_objects["hat_matrices_mv"];
  field<sp_mat> basis_mats_mv = model_objects["basis_matrices_mv"];

  arma::field<cube> V = model_objects["V"];
  arma::field<cube> E = model_objects["E"];
  arma::field<cube> k = model_objects["k"];
  arma::field<cube> eta = model_objects["eta"];
  arma::field<cube> R = model_objects["R"];
  arma::field<cube> R_reg = model_objects["R_reg"];
  arma::field<cube> beta = model_objects["beta"];
  arma::field<cube> beta0field = model_objects["beta0field"];

  //   // Misc parameters
  const int lead_time = model_parameters["lead_time"];
  const std::string loss_function = model_parameters["loss_function"];
  const double loss_parameter = model_parameters["loss_parameter"];
  const bool loss_gradient = model_parameters["loss_gradient"];
  const std::string method = model_parameters["method"];

  const double forget_past_performance = model_parameters["forget_past_performance"];
  const bool allow_quantile_crossing = model_parameters["allow_quantile_crossing"];

  const int start = T - new_y.n_rows;

  cube loss_for(T, D, P, fill::zeros);
  loss_for.rows(0, old_T - 1) =
      Rcpp::as<arma::cube>(object["forecaster_loss"]);
  field<cube> loss_exp(T);
  loss_exp.rows(0, old_T - 1) =
      Rcpp::as<arma::field<cube>>(object["experts_loss"]);

  clock.tock("init");

  clock.tick("core");

  online_learning_core(
      T,
      D,
      P,
      K,
      T_E_Y,
      start,
      lead_time,
      y,
      experts_new,
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
      hat_mats_mv,
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
      prog,
      clock);

  clock.tock("core");

  //   // Update internal objects
  model_objects["V"] = V;
  model_objects["E"] = E;
  model_objects["k"] = k;
  model_objects["eta"] = eta;
  model_objects["R"] = R;
  model_objects["R_reg"] = R_reg;
  model_objects["beta"] = beta;
  model_objects["weights_tmp"] = weights_tmp;
  model_objects["predictions_tmp"] = predictions_tmp;
  model_objects["cum_performance"] = cum_performance;

  //   // Update data
  model_data["experts"] = experts_new;
  model_data["y"] = y;

  Rcpp::NumericMatrix chosen_parameters = Rcpp::wrap(chosen_params);

  object["chosen_parameters"] = chosen_parameters;
  object["opt_index"] = opt_index + 1;
  object["predictions"] = predictions;
  object["weights"] = weights;
  object["past_performance"] = past_performance;
  object["forecaster_loss"] = loss_for;
  object["experts_loss"] = loss_exp;

  return object;
}