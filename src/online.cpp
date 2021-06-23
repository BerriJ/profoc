// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <online.h>
#include <misc.h>
#include <loss.h>
#include <splines2.h>
#include <splines.h>

#include <RcppArmadillo.h>
#include <progress.hpp>

using namespace arma;

//' @template function_online
//'
//' @template param_y
//' @template param_experts
//' @template param_tau
//' @template param_intercept
//' @template param_loss_function
//' @template param_loss_parameter
//' @template param_ex_post_smooth
//' @template param_ex_post_fs
//' @template param_lambda
//' @template param_method
//' @param method_var Allows to calculate slight variations of the BOA
//' algorithm
//' @param forget_regret Share of past regret not to be considered, resp. to be
//' forgotten in every iteration of the algorithm. Defaults to 0.
//' @template param_forget_performance
//' @template param_fixed_share
//' @param gamma Scaling parameter for the learning rate.
//' @template param_ndiff
//' @template param_deg
//' @template param_knot_distance
//' @template param_knot_distance_power
//' @param gradient Determines if a linearized version of the loss is used.
//' @param loss_array User specified loss array. If specified, the loss will not be calculated by profoc.
//' @param regret_array User specified regret array. If specifiec, the regret will not be calculated by profoc.
//' @template param_trace
//' @param init_weights Matrix of dimension Kx1 or KxP used as starting weights. Kx1 represents the constant solution with equal weights over all P whereas specifiying a KxP matrix allows different starting weights for each P.
//' @template param_lead_time
//' @template param_allow_quantile_crossing
//' @template param_soft_threshold
//' @template param_ex_post_soft_threshold
//' @template param_hard_threshold
//' @template param_ex_post_hard_threshold
//' @usage online(y, experts, tau, intercept = FALSE, loss_function = "quantile",
//' loss_parameter = 1, ex_post_smooth = FALSE, ex_post_fs = FALSE,
//' lambda = -Inf, method = "boa", method_var = "A", forget_regret = 0,
//' forget_performance = 0, fixed_share = 0, gamma = 1, ndiff = 1, deg = 3,
//' knot_distance = 0.025, knot_distance_power = 1,
//' gradient = TRUE, loss_array = NULL, regret_array = NULL,
//' trace = TRUE, init_weights = NULL, lead_time = 0, allow_quantile_crossing = FALSE,
//' soft_threshold = -Inf, ex_post_soft_threshold = FALSE, hard_threshold = -Inf,
//' ex_post_hard_threshold = FALSE)
//' @export
// [[Rcpp::export]]
Rcpp::List online(
    mat &y,
    cube &experts,
    Rcpp::NumericVector tau = Rcpp::NumericVector::create(),
    const bool &intercept = false,
    const std::string loss_function = "quantile",
    const double &loss_parameter = 1,
    const bool &ex_post_smooth = false,
    const bool &ex_post_fs = false,
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(-1 / 0),
    const std::string method = "boa",
    const std::string method_var = "A",
    Rcpp::NumericVector forget_regret = Rcpp::NumericVector::create(0),
    const double &forget_performance = 0,
    Rcpp::NumericVector fixed_share = Rcpp::NumericVector::create(0),
    Rcpp::NumericVector gamma = Rcpp::NumericVector::create(1),
    Rcpp::NumericVector ndiff = Rcpp::NumericVector::create(1.5),
    Rcpp::NumericVector deg = Rcpp::NumericVector::create(3),
    Rcpp::NumericVector knot_distance = Rcpp::NumericVector::create(0.025),
    Rcpp::NumericVector knot_distance_power = Rcpp::NumericVector::create(1),
    const bool &gradient = true,
    Rcpp::NumericVector loss_array = Rcpp::NumericVector::create(),
    Rcpp::NumericVector regret_array = Rcpp::NumericVector::create(),
    const bool trace = true,
    Rcpp::Nullable<Rcpp::NumericMatrix> init_weights = R_NilValue,
    const int &lead_time = 0,
    bool allow_quantile_crossing = false,
    Rcpp::NumericVector soft_threshold = Rcpp::NumericVector::create(-1 / 0),
    bool ex_post_soft_threshold = false,
    Rcpp::NumericVector hard_threshold = Rcpp::NumericVector::create(-1 / 0),
    bool ex_post_hard_threshold = false)
{

  if (intercept)
  {
    mat intercept_mat(experts.n_rows, experts.n_cols, fill::ones);
    experts = join_slices(intercept_mat, experts);
  }

  // Indexing Convention -> (T, P, K, X)
  // T number of observations
  // P lengths of probability Grid
  // K number of experts
  // X number of parameter combinations to consider

  // Object Dimensions
  const int T = y.n_rows;
  const int P = experts.n_cols;
  const int K = experts.n_slices;
  const int T_E_Y = experts.n_rows - y.n_rows;

  if (y.n_cols > 1 && !allow_quantile_crossing)
  {
    Rcpp::warning("Warning: allow_quantile_crossing set to true since multivariate prediction target was provided.");
    allow_quantile_crossing = true;
  }

  // Expand y matrix if necessary
  if (y.n_cols == 1)
  {
    y = repmat(y, 1, P);
  }

  // Set default value to tau and / or expand if necessary
  vec tau_vec(tau);
  if (tau_vec.size() == 0)
  {
    tau_vec.resize(P);
    tau_vec = regspace(1, P) / (P + 1);
  }
  else if (tau_vec.size() == 1)
  {
    tau_vec.resize(P);
    tau_vec.fill(tau_vec(0));
  }

  // Init parametergrid
  mat param_grid = get_combinations(lambda, forget_regret);       // Index 0 & 1
  param_grid = get_combinations(param_grid, fixed_share);         // index 2
  param_grid = get_combinations(param_grid, gamma);               // Index 3
  param_grid = get_combinations(param_grid, knot_distance);       // Index 4
  param_grid = get_combinations(param_grid, deg);                 // Index 5
  param_grid = get_combinations(param_grid, ndiff);               // index 6
  param_grid = get_combinations(param_grid, knot_distance_power); // Index 7
  param_grid = get_combinations(param_grid, soft_threshold);      // Index 8
  param_grid = get_combinations(param_grid, hard_threshold);      // Index 9

  const int X = param_grid.n_rows;
  mat chosen_params(T, param_grid.n_cols);
  vec opt_index(T + 1, fill::zeros);
  cube past_performance(T, P, X, fill::zeros);
  vec tmp_performance(X, fill::zeros);
  vec cum_performance(X, fill::zeros);
  Progress prog(T * X + X, trace);

  // Init weight objects

  mat w0;
  // Populate uniform weights if w0 was not specified
  if (init_weights.isNotNull())
  {
    w0 = Rcpp::as<arma::mat>(init_weights);
  }
  else
  {
    w0.resize(K);
    w0.fill(1 / double(K - intercept));
    if (intercept)
      w0.row(0).fill(0);
  }

  // Expand w0 if necessary
  if (w0.n_cols == 1)
  {
    w0 = repmat(w0, 1, P);
  }

  // Truncate from below
  w0.rows(intercept, w0.n_rows - 1) =
      pmax_arma(w0.rows(intercept, w0.n_rows - 1), exp(-350));

  // Normalize weights
  w0 = w0.each_row() / sum(w0);

  // Init object holding temp. weights, resp. ex-ante
  cube w_temp(P, K, X);
  w_temp.each_slice() = w0.t(); // TODO change orientation of w0
  // Init object holding final weights, resp. ex-post
  cube w_post(w_temp);

  // Weights output
  cube weights(T + 1, P, K, fill::zeros);

  cube R(P, K, X, fill::zeros);
  cube R_reg(R);

  cube V(P, K, X, fill::zeros);
  cube E(P, K, X, fill::zeros);
  cube k(P, K, X, fill::zeros);
  cube eta(P, K, X, fill::zeros);
  if (method == "ml_poly")
    eta.fill(exp(350));

  k = k.fill(-699);
  mat experts_mat(P, K);

  cube predictions_post(T, P, X);
  cube predictions_ante(T, P, X);
  mat predictions_temp(P, K);
  mat predictions_final(T + T_E_Y, P, fill::zeros);

  vec lpred(1);
  vec lpred_new(P);
  vec lexp(K);
  mat lexp_new(P, K);
  mat Q(P, K);
  vec r(K);
  vec r_reg(K);
  cube loss_cube(T, P, K, fill::zeros);

  // Init loss array (if existent)
  if (loss_array.size() > 0)
  {
    vec dim = loss_array.attr("dim");
    cube tmp_cube(loss_array.begin(), dim(0), dim(1), dim(2), false);
    loss_cube = tmp_cube;
  }

  cube regret_cube(T, P, K, fill::zeros);
  // Init regret array (if existent)
  if (regret_array.size() > 0)
  {
    vec dim = regret_array.attr("dim");
    cube tmp_cube(regret_array.begin(), dim(0), dim(1), dim(2), false);
    regret_cube = tmp_cube;
  }

  // Smoothing Setup
  field<mat> hat_mats(param_grid.n_rows);
  field<mat> basis_mats(param_grid.n_rows);
  vec spline_basis_x = regspace(1, P) / (P + 1);
  field<mat> Vfield(param_grid.n_rows);
  field<mat> Efield(param_grid.n_rows);
  field<mat> kfield(param_grid.n_rows);
  field<mat> etafield(param_grid.n_rows);
  field<mat> Rfield(param_grid.n_rows);
  field<mat> R_regfield(param_grid.n_rows);
  field<mat> beta(param_grid.n_rows);
  field<mat> w0field(param_grid.n_rows);

  // Init hat matrix field
  for (unsigned int x = 0; x < X; x++)
  {
    // In second step: skip if recalc is not necessary:
    if (x > 0 &&
        param_grid(x, 4) == param_grid(x - 1, 4) &&
        param_grid(x, 5) == param_grid(x - 1, 5) &&
        param_grid(x, 7) == param_grid(x - 1, 7))
    {
      basis_mats(x) = basis_mats(x - 1);
    }
    else
    {

      basis_mats(x) = make_basis_matrix(spline_basis_x,
                                        param_grid(x, 4), // kstep
                                        param_grid(x, 5), // degree
                                        param_grid(x, 7)  // uneven grid
      );
    }

    int L = basis_mats(x).n_cols;
    Vfield(x).zeros(L, K);
    Efield(x).zeros(L, K);
    kfield(x).zeros(L, K);

    mat eta_(L, K, fill::zeros);
    etafield(x) = eta_;
    if (method == "ml_poly")
    {
      eta_.fill(exp(350));
      etafield(x) = eta_;
    }

    R_regfield(x).zeros(L, K);
    Rfield(x).zeros(L, K);

    beta(x) = (w0 * pinv(basis_mats(x)).t()).t();

    w0field(x) = beta(x);

    R_CheckUserInterrupt();
  }

  // Only if smoothing is possible (tau_vec.size > 1)
  if (P > 1)
  {
    // Init hat matrix field
    for (unsigned int x = 0; x < X; x++)
    {
      // In second step: skip if recalc is not necessary:
      if (x > 0 &&
          param_grid(x, 0) == param_grid(x - 1, 0) &&
          param_grid(x, 4) == param_grid(x - 1, 4) &&
          param_grid(x, 5) == param_grid(x - 1, 5) &&
          param_grid(x, 6) == param_grid(x - 1, 6) &&
          param_grid(x, 7) == param_grid(x - 1, 7))
      {
        hat_mats(x) = hat_mats(x - 1);
      }
      else
      {

        hat_mats(x) = make_hat_matrix(spline_basis_x,
                                      param_grid(x, 4), // kstep
                                      param_grid(x, 0), // lambda
                                      param_grid(x, 6), // differences
                                      param_grid(x, 5), // degree
                                      param_grid(x, 7)  // uneven grid
        );
      }
      R_CheckUserInterrupt();
      prog.increment(); // Update progress
    }
  }

  // Predictions at t < lead_time using initial weights
  for (unsigned int t = 0; t < lead_time; t++)
  {
    // Save final weights w_post
    weights.row(t) = w_post.slice(opt_index(t));

    // Store expert predictions temporarily
    experts_mat = experts.row(t);

    // Forecasters prediction(ex-ante)
    predictions_temp = sum(w_temp.each_slice() % experts_mat, 1);
    predictions_ante.row(t) = predictions_temp;

    // Forecasters prediction (ex-post)
    // Note that w_post != w_temp in ex post setting and w_post = w_temp otherwise
    predictions_temp = sum(w_post.each_slice() % experts_mat, 1);
    predictions_post.row(t) = predictions_temp;

    // Final prediction
    predictions_final.row(t) =
        vectorise(predictions_post(span(t), span::all, span(opt_index(t)))).t();

    past_performance.row(t).fill(datum::nan);
  }

  for (unsigned int t = (0 + lead_time); t < T; t++)
  {

    // Save final weights w_post
    weights.row(t) = w_post.slice(opt_index(t));

    // Store expert predictions temporarily
    experts_mat = experts.row(t);

    // Forecasters prediction(ex-ante)
    predictions_temp = sum(w_temp.each_slice() % experts_mat, 1);
    predictions_ante.row(t) = predictions_temp;

    // Forecasters prediction (ex-post)
    // Note that w_post != w_temp in ex post setting and w_post = w_temp otherwise
    predictions_temp = sum(w_post.each_slice() % experts_mat, 1);
    predictions_post.row(t) = predictions_temp;

    // Final prediction
    predictions_final.row(t) =
        vectorise(predictions_post(span(t), span::all, span(opt_index(t)))).t();

    for (unsigned int x = 0; x < X; x++)
    {

      // NEW -------------------------------------------------------------------
      for (unsigned int p = 0; p < P; p++)
      {

        // Evaluate the ex-post predictive performance
        past_performance(t, p, x) = loss(y(t, p),
                                         predictions_post(t - lead_time, p, x),
                                         9999,           // where evaluate gradient
                                         loss_function,  // method
                                         tau_vec(p),     // tau_vec
                                         loss_parameter, // alpha
                                         false);

        for (unsigned int k = 0; k < K; k++)
        {
          lexp_new(p, k) = loss(y(t, p),
                                experts(t, p, k),
                                predictions_ante(t - lead_time, p, x), // where evaluate gradient
                                loss_function,                         // method
                                tau_vec(p),                            // tau_vec
                                loss_parameter,                        // alpha
                                gradient);
        }

        if (loss_array.size() != 0)
        {
          lexp_new.row(p) = vectorise(loss_cube.tube(t, p));
        }

        lpred_new(p) = loss(y(t, p),
                            predictions_ante(t - lead_time, p, x),
                            predictions_ante(t - lead_time, p, x), // where to evaluate gradient
                            loss_function,                         // method
                            tau_vec(p),                            // tau_vec
                            loss_parameter,                        // alpha
                            gradient);

        Q.row(p) = lpred_new(p) - lexp_new.row(p);
      }

      mat Q_regret;

      if (regret_array.size() == 0)
      {
        Q_regret = Q.t() * basis_mats(x);
      }
      else
      {
        Q_regret = regret_cube.row(t);
      }

      for (unsigned int l = 0; l < Q_regret.n_cols; l++)
      {

        r = Q_regret.col(l);

        if (method == "ewa")
        {
          // Update the cumulative regret used by eta
          Rfield(x).row(l) *= (1 - param_grid(x, 1));
          Rfield(x).row(l) += r.t();

          beta(x).row(l) = exp(param_grid(x, 3) * Rfield(x).row(l));
          beta(x).row(l) /= accu(beta(x).row(l));
        }
        else if (method == "ml_poly")
        {
          // Update the cumulative regret used by ML_Poly
          Rfield(x).row(l) *= (1 - param_grid(x, 1));
          Rfield(x).row(l) += r.t();

          //   // Update the learning rate
          etafield(x).row(l) = 1 / (1 / etafield(x).row(l) + square(r.t()));

          // param_grid(x, 3) = gamma
          beta(x).row(l) =
              param_grid(x, 3) * etafield(x).row(l) % pmax_arma(Rfield(x).row(l), exp(-700));
          beta(x).row(l) /= accu(beta(x).row(l));
        }
        else if (method == "boa")
        {

          Vfield(x).row(l) *= (1 - param_grid(x, 1));
          Vfield(x).row(l) += square(r.t());

          Efield(x).row(l) =
              max(Efield(x).row(l) * (1 - param_grid(x, 1)), abs(r.t()));

          etafield(x).row(l) =
              pmin_arma(
                  min(1 / (2 * Efield(x).row(l)),
                      sqrt(log(K) / Vfield(x).row(l))),
                  exp(350));

          r_reg = r - etafield(x).row(l).t() % square(r);

          if (method_var.find('A') != std::string::npos)
          {
            R_regfield(x).row(l) *= (1 - param_grid(x, 1));
            R_regfield(x).row(l) +=
                0.5 * (r_reg.t() + conv_to<colvec>::from(etafield(x).row(l) % r.t() > 0.5).t() % (2 * Efield(x).row(l)));
          }
          else
          {
            // Gaillard
            R_regfield(x).row(l) *= (1 - param_grid(x, 1));
            R_regfield(x).row(l) += r_reg.t();
          }

          if (method_var.find('C') != std::string::npos)
          {
            // Wintenberger
            beta(x).row(l) =
                param_grid(x, 3) * etafield(x).row(l) % exp(param_grid(x, 3) * etafield(x).row(l) % R_regfield(x).row(l)) % w0field(x).row(l);

            beta(x).row(l) /=
                mean(param_grid(x, 3) * etafield(x).row(l) % exp(param_grid(x, 3) * etafield(x).row(l) % R_regfield(x).row(l)));
          }
          else
          {
            // Gaillard
            beta(x).row(l) =
                exp(param_grid(x, 3) * etafield(x).row(l) % R_regfield(x).row(l));
            beta(x).row(l) = pmin_arma(pmax_arma(beta(x).row(l), exp(-700)), exp(700));

            beta(x).row(l) /= accu(beta(x).row(l));
          }
        }
        else
        {
          Rcpp::stop("Choose 'boa', 'ml_poly' or 'ewa' as method.");
        }
      }

      w_temp.slice(x) = basis_mats(x) * beta(x);

      w_post.slice(x) = w_temp.slice(x);

      // Apply thresholds

      if (!ex_post_soft_threshold)
      {
        for (double &e : w_temp.slice(x))
        {
          threshold_soft(e, param_grid(x, 8));
        }
      }

      for (double &e : w_post.slice(x))
      {
        threshold_soft(e, param_grid(x, 8));
      }

      if (!ex_post_hard_threshold)
      {
        for (double &e : w_temp.slice(x))
        {
          threshold_hard(e, param_grid(x, 9));
        }
      }

      for (double &e : w_post.slice(x))
      {
        threshold_hard(e, param_grid(x, 9));
      }

      // Smoothing
      if (param_grid(x, 0) != -datum::inf)
      {
        for (unsigned int k = 0; k < K; k++)
        {
          if (!ex_post_smooth)
          {
            w_temp(span::all, span(k), span(x)) =
                hat_mats(x) * vectorise(w_temp(span::all, span(k), span(x)));
          }

          w_post(span::all, span(k), span(x)) = hat_mats(x) * vectorise(w_temp(span::all, span(k), span(x)));
        }
      }

      for (unsigned int p = 0; p < P; p++)
      {
        //Add fixed_share

        if (!ex_post_fs)
        {
          w_temp(span(p), span::all, span(x)) =
              (1 - param_grid(x, 2)) * w_temp(span(p), span::all, span(x)) +
              param_grid(x, 2) / K;
        }

        w_post(span(p), span::all, span(x)) =
            (1 - param_grid(x, 2)) * w_post(span(p), span::all, span(x)) +
            (param_grid(x, 2) / K);

        // Enshure that constraints hold

        // Positivity
        w_temp(span(p), span::all, span(x)) =
            pmax_arma(w_temp(span(p), span::all, span(x)), exp(-700));
        w_post(span(p), span::all, span(x)) =
            pmax_arma(w_post(span(p), span::all, span(x)), exp(-700));

        // Affinity
        w_post(span(p), span::all, span(x)) /=
            accu(w_post(span(p), span::all, span(x)));
        w_temp(span(p), span::all, span(x)) /=
            accu(w_temp(span(p), span::all, span(x)));
      }

      tmp_performance(x) = accu(past_performance(span(t), span::all, span(x)));
      prog.increment(); // Update progress
      R_CheckUserInterrupt();
    }

    // Sum past_performance in each slice
    cum_performance = (1 - forget_performance) * cum_performance + tmp_performance;
    opt_index(t + 1) = cum_performance.index_min();
    chosen_params.row(t) = param_grid.row(opt_index(t + 1));
  }

  // Save Weights and Prediction
  weights.row(T) = w_post.slice(opt_index(T));

  // Predict residual expert forecasts if any are available
  for (unsigned int t = T; t < T + T_E_Y; t++)
  {
    experts_mat = experts.row(t);
    predictions_temp = sum(w_post.slice(opt_index(T)) % experts_mat, 1);
    predictions_final.row(t) = vectorise(predictions_temp).t();
  }

  // Sort predictions if quantile_crossing is prohibited
  if (!allow_quantile_crossing)
  {
    predictions_final = arma::sort(predictions_final, "ascend", 1);
  }

  // Save losses suffered by forecaster and experts
  mat loss_for(T, P, fill::zeros);
  cube loss_exp(T, P, K, fill::zeros);

  for (unsigned int t = 0; t < T; t++)
  {
    for (unsigned int p = 0; p < P; p++)
    {
      for (unsigned int k = 0; k < K; k++)
      {
        loss_exp(t, p, k) =
            loss(y(t, p),
                 experts(t, p, k),
                 9999,           // where to evaluate the gradient
                 loss_function,  // method
                 tau_vec(p),     // tau_vec
                 loss_parameter, // alpha
                 false);         // gradient
      }
      loss_for(t, p) = loss(y(t, p),
                            predictions_final(t, p),
                            9999,           // where to evaluate the gradient
                            loss_function,  // method
                            tau_vec(p),     // tau_vec
                            loss_parameter, // alpha
                            false);         // gradient;
    }
  }

  // 1-Indexing for R-Output
  opt_index = opt_index + 1;

  // Set unused values to NA
  if (lead_time > 0)
  {
    chosen_params.rows(0, lead_time - 1).fill(datum::nan);
  }

  Rcpp::DataFrame chosen_parameters = Rcpp::DataFrame::create(
      Rcpp::Named("lambda") = chosen_params.col(0),
      Rcpp::Named("forget_regret") = chosen_params.col(1),
      Rcpp::Named("fixed_share") = chosen_params.col(2),
      Rcpp::Named("gamma") = chosen_params.col(3),
      Rcpp::Named("knot_distance") = chosen_params.col(4),
      Rcpp::Named("deg") = chosen_params.col(5),
      Rcpp::Named("diff") = chosen_params.col(6),
      Rcpp::Named("knot_distance_power") = chosen_params.col(7),
      Rcpp::Named("threshold_soft") = chosen_params.col(8),
      Rcpp::Named("threshold_hard") = chosen_params.col(9));

  Rcpp::DataFrame parametergrid = Rcpp::DataFrame::create(
      Rcpp::Named("lambda") = param_grid.col(0),
      Rcpp::Named("forget_regret") = param_grid.col(1),
      Rcpp::Named("fixed_share") = param_grid.col(2),
      Rcpp::Named("gamma") = param_grid.col(3),
      Rcpp::Named("knot_distance") = param_grid.col(4),
      Rcpp::Named("deg") = param_grid.col(5),
      Rcpp::Named("diff") = param_grid.col(6),
      Rcpp::Named("knot_distance_power") = param_grid.col(7),
      Rcpp::Named("threshold_soft") = param_grid.col(8),
      Rcpp::Named("threshold_hard") = param_grid.col(9));

  return Rcpp::List::create(
      Rcpp::Named("predictions") = predictions_final,
      Rcpp::Named("weights") = weights,
      Rcpp::Named("forecaster_loss") = loss_for,
      Rcpp::Named("experts_loss") = loss_exp,
      Rcpp::Named("past_perf_wrt_params") = past_performance,
      Rcpp::Named("chosen_parameters") = chosen_parameters,
      Rcpp::Named("parametergrid") = parametergrid,
      Rcpp::Named("opt_index") = opt_index);
}
