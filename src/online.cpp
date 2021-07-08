// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <misc.h>
#include <loss.h>
#include <splines2.h>
#include <splines.h>

#include <RcppArmadillo.h>
#include <progress.hpp>

using namespace arma;

// TODO Move core objects into a struct

void online_learning_core(
    const int &T,
    const int &P,
    const int &K,
    const int &T_E_Y,
    const int &start,
    const int &lead_time,
    const mat &y,
    const cube &experts,
    const vec &tau_vec,
    const std::string loss_function,
    const std::string method,
    const std::string method_var,
    const double &loss_parameter,
    const bool &loss_gradient,
    cube &weights,
    cube &weights_tmp,
    cube &predictions_tmp,
    mat &predictions,
    cube &past_performance,
    vec &opt_index,
    mat &param_grid,
    mat &chosen_params,
    field<mat> &R,
    field<mat> &R_reg,
    field<mat> &V,
    field<mat> &E,
    field<mat> &eta,
    field<mat> &hat_mats,
    field<sp_mat> &basis_mats,
    field<mat> &w0field,
    field<mat> &beta,
    vec &cum_performance,
    const double &forget_past_performance,
    vec &tmp_performance,
    const bool &allow_quantile_crossing,
    mat &loss_for,
    cube &loss_exp,
    cube &loss_cube,
    cube &regret_cube,
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
    predictions_tmp.row(t) = predictions_temp;

    // Final prediction
    predictions.row(t) =
        vectorise(predictions_tmp(span(t), span::all, span(opt_index(t)))).t();

    for (unsigned int x = 0; x < param_grid.n_rows; x++)
    {

      mat lexp(P, K); // Experts loss
      vec lfor(P);    // Forecasters loss

      for (unsigned int p = 0; p < P; p++)
      {

        // Evaluate the ex-post predictive performance
        past_performance(t, p, x) = loss(y(t, p),
                                         predictions_tmp(t - lead_time, p, x),
                                         9999,           // where evaluate loss_gradient
                                         loss_function,  // method
                                         tau_vec(p),     // tau_vec
                                         loss_parameter, // alpha
                                         false);

        for (unsigned int k = 0; k < K; k++)
        {
          lexp(p, k) = loss(y(t, p),
                            experts(t, p, k),
                            predictions_tmp(t - lead_time, p, x), // where evaluate loss_gradient
                            loss_function,                        // method
                            tau_vec(p),                           // tau_vec
                            loss_parameter,                       // alpha
                            loss_gradient);
        }

        if (loss_cube.n_elem != 0)
        {
          lexp.row(p) = vectorise(loss_cube.tube(t, p)).t();
        }

        lfor(p) = loss(y(t, p),
                       predictions_tmp(t - lead_time, p, x),
                       predictions_tmp(t - lead_time, p, x), // where to evaluate loss_gradient
                       loss_function,                        // method
                       tau_vec(p),                           // tau_vec
                       loss_parameter,                       // alpha
                       loss_gradient);
      }

      mat Q_regret;
      if (regret_cube.n_elem == 0)
      {
        Q_regret = (lfor - lexp.each_col()).t();
        Q_regret *= double(basis_mats(x).n_cols) / double(P);
        Q_regret *= basis_mats(x);
      }
      else
      {
        Q_regret = regret_cube.row(t);
        Q_regret = Q_regret.t();
        Q_regret *= double(basis_mats(x).n_cols) / double(P);
        Q_regret *= basis_mats(x);
      }

      for (unsigned int l = 0; l < Q_regret.n_cols; l++)
      {

        vec r = Q_regret.col(l);

        if (method == "ewa")
        {
          // Update the cumulative regret used by eta
          R(x).row(l) *= (1 - param_grid(x, 3));
          R(x).row(l) += r.t();

          beta(x).row(l) = exp(param_grid(x, 12) * R(x).row(l));
          beta(x).row(l) /= accu(beta(x).row(l));
        }
        else if (method == "ml_poly")
        {
          // Update the cumulative regret used by ML_Poly
          R(x).row(l) *= (1 - param_grid(x, 3));
          R(x).row(l) += r.t();

          //   // Update the learning rate
          eta(x).row(l) = 1 / (1 / eta(x).row(l) + square(r.t()));

          // param_grid(x, 12) = gamma
          beta(x).row(l) =
              param_grid(x, 12) * eta(x).row(l) % pmax_arma(R(x).row(l), exp(-700));
          beta(x).row(l) /= accu(beta(x).row(l));
        }
        else if (method == "boa")
        {

          V(x).row(l) *= (1 - param_grid(x, 3));
          V(x).row(l) += square(r.t());

          E(x).row(l) =
              max(E(x).row(l) * (1 - param_grid(x, 3)), abs(r.t()));

          eta(x).row(l) =
              pmin_arma(
                  min(1 / (2 * E(x).row(l)),
                      sqrt(log(K) / V(x).row(l))),
                  exp(350));

          vec r_reg = r - eta(x).row(l).t() % square(r);

          if (method_var.find('A') != std::string::npos)
          {
            R_reg(x).row(l) *= (1 - param_grid(x, 3));
            R_reg(x).row(l) +=
                0.5 * (r_reg.t() + conv_to<colvec>::from(eta(x).row(l) % r.t() > 0.5).t() % (2 * E(x).row(l)));
          }
          else
          {
            // Gaillard
            R_reg(x).row(l) *= (1 - param_grid(x, 3));
            R_reg(x).row(l) += r_reg.t();
          }

          if (method_var.find('C') != std::string::npos)
          {
            // Wintenberger
            beta(x).row(l) =
                param_grid(x, 12) * eta(x).row(l) % exp(param_grid(x, 12) * eta(x).row(l) % R_reg(x).row(l)) % w0field(x).row(l);

            beta(x).row(l) /=
                mean(param_grid(x, 12) * eta(x).row(l) % exp(param_grid(x, 12) * eta(x).row(l) % R_reg(x).row(l)));
          }
          else
          {
            // Gaillard
            beta(x).row(l) =
                exp(param_grid(x, 12) * eta(x).row(l) % R_reg(x).row(l));
            beta(x).row(l) = pmin_arma(pmax_arma(beta(x).row(l), exp(-700)), exp(700));

            beta(x).row(l) /= accu(beta(x).row(l));
          }
        }
        else
        {
          Rcpp::stop("Choose 'boa', 'ml_poly' or 'ewa' as method.");
        }

        // Apply thresholds
        for (double &e : beta(x).row(l))
        {
          threshold_soft(e, param_grid(x, 4));
        }

        for (double &e : beta(x).row(l))
        {
          threshold_hard(e, param_grid(x, 5));
        }

        //Add fixed_share
        beta(x).row(l) =
            (1 - param_grid(x, 6)) * beta(x).row(l) +
            (param_grid(x, 6) / K);
      }

      weights_tmp.slice(x) = basis_mats(x) * beta(x);

      // Smoothing
      if (param_grid(x, 7) != -datum::inf)
      {
        for (unsigned int k = 0; k < K; k++)
        {
          weights_tmp(span::all, span(k), span(x)) = hat_mats(x) * vectorise(weights_tmp(span::all, span(k), span(x)));
        }
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
    predictions.row(t) = vectorise(predictions_temp).t();
  }

  // Sort predictions if quantile_crossing is prohibited
  if (!allow_quantile_crossing)
  {
    predictions = arma::sort(predictions, "ascend", 1);
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
                 tau_vec(p),     // tau_vec
                 loss_parameter, // alpha
                 false);         // loss_gradient
      }
      loss_for(t, p) = loss(y(t, p),
                            predictions(t, p),
                            9999,           // where to evaluate the loss_gradient
                            loss_function,  // method
                            tau_vec(p),     // tau_vec
                            loss_parameter, // alpha
                            false);         // loss_gradient;
    }
  }

  return;
}

//' @template function_online
//'
//' @template param_y
//' @template param_experts
//' @template param_tau
//'
//' @template param_lead_time
//'
//' @template param_loss_function
//' @template param_loss_parameter
//' @param loss_gradient Determines if a linearized version of the loss is used.
//'
//' @template param_method
//' @param method_var Allows to calculate slight variations of the BOA
//' algorithm
//'
//' @template param_basis_knot_distance_online
//' @template param_basis_knot_distance_power
//' @template param_basis_deg_online
//'
//' @param forget_regret Share of past regret not to be considered, resp. to be
//' forgotten in every iteration of the algorithm. Defaults to 0.
//'
//' @template param_soft_threshold
//' @template param_hard_threshold
//'
//' @template param_fixed_share
//'
//' @template param_smooth_lambda
//' @template param_smooth_knot_distance
//' @template param_smooth_knot_distance_power
//' @template param_smooth_deg
//' @template param_smooth_ndiff
//'
//' @param gamma Scaling parameter for the learning rate.
//'
//' @template param_parametergrid_max_combinations
//' @template param_parametergrid
//' @template param_forget_past_performance
//'
//' @template param_allow_quantile_crossing
//'
//' @param init_weights Matrix of dimension 1xK or PxK used as starting weights. 1xK represents the constant solution with equal weights over all P whereas specifiying a PxK matrix allows different starting weights for each P.
//' @param loss_array User specified loss array. If specified, the loss will not be calculated by profoc.
//' @param regret_array User specified regret array. If specifiec, the regret will not be calculated by profoc.
//' @template param_trace
//'
//' @usage online(
//' y,
//' experts,
//' tau,
//' lead_time = 0,
//' loss_function = "quantile",
//' loss_parameter = 1,
//' loss_gradient = TRUE,
//' method = "boa",
//' method_var = "A",
//' basis_knot_distance = c(2^seq(log(1/(length(tau)+1),2)-1, -1, length=5),1),
//' basis_knot_distance_power = 1,
//' basis_deg = 3,
//' forget_regret = 0,
//' soft_threshold = -Inf,
//' hard_threshold = -Inf,
//' fixed_share = 0,
//' smooth_lambda = -Inf,
//' smooth_knot_distance = c(2^seq(log(1/(length(tau)+1),2)-1, -1, length=5),1),
//' smooth_knot_distance_power = 1,
//' smooth_deg = 3,
//' smooth_ndiff = 1.5,
//' gamma = 1,
//' parametergrid_max_combinations = 100,
//' parametergrid = NULL,
//' forget_past_performance = 0,
//' allow_quantile_crossing = FALSE,
//' init_weights = NULL,
//' loss_array = NULL,
//' regret_array = NULL,
//' trace = TRUE
//' )
//' @export
// [[Rcpp::export]]
Rcpp::List online(
    mat &y,
    cube &experts,
    Rcpp::NumericVector tau = Rcpp::NumericVector::create(),
    const int &lead_time = 0,
    const std::string loss_function = "quantile",
    const double &loss_parameter = 1,
    const bool &loss_gradient = true,
    const std::string method = "boa",
    const std::string method_var = "A",
    Rcpp::NumericVector basis_knot_distance = Rcpp::NumericVector::create(),
    Rcpp::NumericVector basis_knot_distance_power = Rcpp::NumericVector::create(1),
    Rcpp::NumericVector basis_deg = Rcpp::NumericVector::create(3),
    Rcpp::NumericVector forget_regret = Rcpp::NumericVector::create(0),
    Rcpp::NumericVector soft_threshold = Rcpp::NumericVector::create(-1 / 0),
    Rcpp::NumericVector hard_threshold = Rcpp::NumericVector::create(-1 / 0),
    Rcpp::NumericVector fixed_share = Rcpp::NumericVector::create(0),
    Rcpp::NumericVector smooth_lambda = Rcpp::NumericVector::create(-1 / 0),
    Rcpp::NumericVector smooth_knot_distance = Rcpp::NumericVector::create(),
    Rcpp::NumericVector smooth_knot_distance_power = Rcpp::NumericVector::create(),
    Rcpp::NumericVector smooth_deg = Rcpp::NumericVector::create(),
    Rcpp::NumericVector smooth_ndiff = Rcpp::NumericVector::create(1.5),
    Rcpp::NumericVector gamma = Rcpp::NumericVector::create(1),
    const int parametergrid_max_combinations = 100,
    Rcpp::Nullable<Rcpp::NumericMatrix> parametergrid = R_NilValue,
    const double &forget_past_performance = 0,
    bool allow_quantile_crossing = false,
    Rcpp::Nullable<Rcpp::NumericMatrix> init_weights = R_NilValue,
    Rcpp::NumericVector loss_array = Rcpp::NumericVector::create(),
    Rcpp::NumericVector regret_array = Rcpp::NumericVector::create(),
    const bool trace = true)
{

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

  if (T_E_Y < 0)
    Rcpp::stop("Number of provided expert predictions has to match or exceed observations.");

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

  vec basis_knot_distance_vec = basis_knot_distance;
  if (basis_knot_distance.size() == 0)
  {
    vec tmp(6, fill::zeros);
    tmp.subvec(0, 4) = arma::linspace(log2(1 / (double(P) + 1)) - 1, -1, 5);
    for (double &e : tmp)
    {
      e = pow(2, e);
    }
    basis_knot_distance_vec = tmp;
  }

  bool inh_deg = false;
  if (smooth_deg.size() == 0)
    inh_deg = true;

  bool inh_kstep = false;
  if (smooth_knot_distance.size() == 0)
    inh_kstep = true;

  bool inh_kstep_p = false;
  if (smooth_knot_distance_power.size() == 0)
    inh_kstep_p = true;

  // Init parametergrid
  mat param_grid;

  if (parametergrid.isNotNull())
  {
    param_grid = Rcpp::as<arma::mat>(parametergrid);
    if (param_grid.n_cols != 13)
      Rcpp::stop("Please provide a parametergrid with 13 columns.");
  }
  else
  {
    param_grid =
        get_combinations(basis_knot_distance_vec,                                          // Index 0
                         basis_knot_distance_power);                                       // Index 1
    param_grid = get_combinations(param_grid, basis_deg);                                  // index 2
    param_grid = get_combinations(param_grid, forget_regret);                              // index 3
    param_grid = get_combinations(param_grid, soft_threshold);                             // index 4
    param_grid = get_combinations(param_grid, hard_threshold);                             // index 5
    param_grid = get_combinations(param_grid, fixed_share);                                // index 6
    param_grid = get_combinations(param_grid, smooth_lambda);                              // index 7
    param_grid = get_combinations(param_grid, smooth_knot_distance, inh_kstep, 0);         // Index 8
    param_grid = get_combinations(param_grid, smooth_knot_distance_power, inh_kstep_p, 1); // Index 9
    param_grid = get_combinations(param_grid, smooth_deg, inh_deg, 2);                     // Index 10
    param_grid = get_combinations(param_grid, smooth_ndiff);                               // Index 11
    param_grid = get_combinations(param_grid, gamma);                                      // Index 12
  }

  if (param_grid.n_rows > parametergrid_max_combinations)
  {
    Rcpp::warning("Warning: Too many parameter combinations possible. %m combinations were randomly sampled. Results may depend on sampling.", parametergrid_max_combinations);
    uvec tmp_index = randperm(param_grid.n_rows, parametergrid_max_combinations);
    tmp_index = sort(tmp_index);
    param_grid = param_grid.rows(tmp_index);
  }

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
    if ((w0.n_rows != 1 && w0.n_rows != P) || w0.n_cols != K)
      Rcpp::stop("Either a 1xK or PxK matrix of initial weights must be supplied.");
  }
  else
  {
    w0.set_size(1, K);
    w0.fill(1 / double(K));
  }

  // Expand w0 if necessary
  if (w0.n_rows == 1)
  {
    w0 = repmat(w0, P, 1);
  }

  // Truncate from below
  w0 = pmax_arma(w0, exp(-350));

  // Normalize weights
  w0.each_col() /= sum(w0, 1);

  // Init object holding weights
  cube weights_tmp(P, K, X);
  weights_tmp.each_slice() = w0;

  cube predictions_tmp(T, P, X);

  // Output Objects
  mat predictions(T + T_E_Y, P, fill::zeros);
  cube weights(T + 1, P, K, fill::zeros);

  // Init loss array (if existent)
  cube loss_cube;
  if (loss_array.size() > 0)
  {
    vec dim = loss_array.attr("dim");
    cube tmp_cube(loss_array.begin(), dim(0), dim(1), dim(2), false);
    loss_cube = tmp_cube;
  }
  // Init regret array (if existent)
  cube regret_cube;
  if (regret_array.size() > 0)
  {
    vec dim = regret_array.attr("dim");
    cube tmp_cube(regret_array.begin(), dim(0), dim(1), dim(2), false);
    regret_cube = tmp_cube;
  }

  // Learning parameters
  field<mat> hat_mats(X);
  field<sp_mat> basis_mats(X);
  field<mat> V(X);
  field<mat> E(X);
  field<mat> k(X);
  field<mat> eta(X);
  field<mat> R(X);
  field<mat> R_reg(X);
  field<mat> beta(X);
  field<mat> w0field(X);

  vec spline_basis_x = regspace(1, P) / (P + 1);

  // Init hat matrix field
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

      basis_mats(x) = make_basis_matrix(spline_basis_x,
                                        param_grid(x, 0),  // kstep
                                        param_grid(x, 2),  // degree
                                        param_grid(x, 1)); // uneven grid
    }

    int L = basis_mats(x).n_cols;
    V(x).zeros(L, K);
    E(x).zeros(L, K);
    k(x).zeros(L, K);

    mat eta_(L, K, fill::zeros);
    eta(x) = eta_;
    if (method == "ml_poly")
    {
      eta_.fill(exp(350));
      eta(x) = eta_;
    }

    R_reg(x).zeros(L, K);
    R(x).zeros(L, K);

    beta(x) = (w0.t() * pinv(mat(basis_mats(x))).t()).t();

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

        hat_mats(x) = make_hat_matrix(spline_basis_x,
                                      param_grid(x, 8),  // kstep
                                      param_grid(x, 7),  // lambda
                                      param_grid(x, 11), // differences
                                      param_grid(x, 10), // degree
                                      param_grid(x, 9)   // uneven grid
        );
      }
      R_CheckUserInterrupt();
      prog.increment(); // Update progress
    }
  }

  if (T <= lead_time)
    Rcpp::stop("Number of expert predictions need to exceed lead_time.");

  // Predictions at t < lead_time using initial weights
  for (unsigned int t = 0; t < lead_time; t++)
  {
    // Save final weights weights_tmp
    weights.row(t) = weights_tmp.slice(opt_index(t));

    // Store expert predictions temporarily
    mat experts_mat = experts.row(t);

    // Forecasters prediction
    mat predictions_temp = sum(weights_tmp.each_slice() % experts_mat, 1);
    predictions_tmp.row(t) = predictions_temp;

    // Final prediction
    predictions.row(t) =
        vectorise(predictions_tmp(span(t), span::all, span(opt_index(t)))).t();

    past_performance.row(t).fill(datum::nan);
  }

  mat loss_for(T, P, fill::zeros);
  cube loss_exp(T, P, K, fill::zeros);

  const int start = lead_time;

  online_learning_core(
      T,
      P,
      K,
      T_E_Y,
      start,
      lead_time,
      y,
      experts,
      tau_vec,
      loss_function,
      method,
      method_var,
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
      w0field,
      beta,
      cum_performance,
      forget_past_performance,
      tmp_performance,
      allow_quantile_crossing,
      loss_for,
      loss_exp,
      loss_cube,
      regret_cube,
      prog);

  // 1-Indexing for R-Output
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
                                    "smooth_lambda",
                                    "smooth_knot_distance",
                                    "smooth_knot_distance_power",
                                    "smooth_deg",
                                    "smooth_diff",
                                    "gamma");

  Rcpp::NumericMatrix parametergrid_out = Rcpp::wrap(param_grid);
  Rcpp::colnames(parametergrid_out) = param_names;

  Rcpp::NumericMatrix chosen_parameters = Rcpp::wrap(chosen_params);
  Rcpp::colnames(chosen_parameters) = param_names;

  Rcpp::List model_data = Rcpp::List::create(
      Rcpp::Named("y") = y,
      Rcpp::Named("experts") = experts,
      Rcpp::Named("tau") = tau_vec);

  Rcpp::List model_parameters = Rcpp::List::create(
      Rcpp::Named("lead_time") = lead_time,
      Rcpp::Named("loss_function") = loss_function,
      Rcpp::Named("loss_parameter") = loss_parameter,
      Rcpp::Named("loss_gradient") = loss_gradient,
      Rcpp::Named("method") = method,
      Rcpp::Named("method_var") = method_var,
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
      Rcpp::Named("w0field") = w0field);

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

  out.attr("class") = "online";

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
    Rcpp::NumericVector new_experts = Rcpp::NumericVector::create())
{

  cube experts_new;

  if (new_experts.size() > 0)
  {
    vec dim = new_experts.attr("dim");
    cube tmp_cube(new_experts.begin(), dim(0), dim(1), dim(2), false);
    experts_new = tmp_cube;
  }

  // This creates a reference to the specification, not a copy
  Rcpp::List specification = object["specification"];

  Rcpp::List model_parameters = specification["parameters"];
  Rcpp::List model_data = specification["data"];
  Rcpp::List model_objects = specification["objects"];

  // Data
  cube experts = model_data["experts"];
  int P = experts.n_cols;

  if (new_experts.size() > 0)
    experts.insert_rows(experts.n_rows, experts_new);

  mat y = model_data["y"];
  // Expand y matrix if necessary
  if (new_y.n_cols == 1)
  {
    new_y = arma::repmat(new_y, 1, P);
  }
  y.insert_rows(y.n_rows, new_y);

  Progress prog(999, false);

  // Object Dimensions
  const int T = y.n_rows;
  const int K = experts.n_slices;
  const int T_E_Y = experts.n_rows - y.n_rows;

  if (T_E_Y < 0)
    Rcpp::stop("Number of provided expert predictions has to match or exceed observations.");

  vec tau_vec = model_data["tau"];
  mat param_grid = object["parametergrid"];
  const int X = param_grid.n_rows;
  mat chosen_params = object["chosen_parameters"];
  chosen_params.resize(T, param_grid.n_cols);
  vec opt_index = object["opt_index"];
  // Zero indexing in C++
  opt_index -= 1;
  opt_index.resize(T + 1);

  cube past_performance = object["past_performance"];
  past_performance.resize(T, P, X);

  vec tmp_performance(X, fill::zeros);
  vec cum_performance = model_objects["cum_performance"];

  cube weights_tmp = model_objects["weights_tmp"];

  cube predictions_tmp = model_objects["predictions_tmp"];
  predictions_tmp.resize(T, P, X);

  cube loss_cube;
  cube regret_cube;

  // Output Objects
  mat predictions = object["predictions"];
  predictions.resize(T + T_E_Y, P);
  cube weights = object["weights"];
  weights.resize(T + 1, P, K);

  // TODO Add loss_cube and regret_array functionality

  field<mat> hat_mats = model_objects["hat_matrices"];
  field<sp_mat> basis_mats = object["basis_matrices"];
  field<mat> V = model_objects["V"];
  field<mat> E = model_objects["E"];
  field<mat> k = model_objects["k"];
  field<mat> eta = model_objects["eta"];
  field<mat> R = model_objects["R"];
  field<mat> R_reg = model_objects["R_reg"];
  field<mat> beta = model_objects["beta"];
  field<mat> w0field = model_objects["w0field"];

  // Misc parameters
  const int lead_time = model_parameters["lead_time"];
  const std::string loss_function = model_parameters["loss_function"];
  const double loss_parameter = model_parameters["loss_parameter"];
  const bool loss_gradient = model_parameters["loss_gradient"];
  const std::string method = model_parameters["method"];
  const std::string method_var = model_parameters["method_var"];

  const double forget_past_performance = model_parameters["forget_past_performance"];
  const bool allow_quantile_crossing = model_parameters["allow_quantile_crossing"];

  const int start = T - new_y.n_rows;

  mat loss_for(T, P, fill::zeros);
  cube loss_exp(T, P, K, fill::zeros);

  online_learning_core(
      T,
      P,
      K,
      T_E_Y,
      start,
      lead_time,
      y,
      experts,
      tau_vec,
      loss_function,
      method,
      method_var,
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
      w0field,
      beta,
      cum_performance,
      forget_past_performance,
      tmp_performance,
      allow_quantile_crossing,
      loss_for,
      loss_exp,
      loss_cube,
      regret_cube,
      prog);

  // Update internal objects
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

  // Update data
  model_data["experts"] = experts;
  model_data["y"] = y;

  // Update output
  Rcpp::CharacterVector param_names =
      Rcpp::CharacterVector::create("basis_knot_distance",
                                    "basis_knot_distance_power",
                                    "basis_deg",
                                    "forget_regret",
                                    "threshold_soft",
                                    "threshold_hard",
                                    "fixed_share",
                                    "smooth_lambda",
                                    "smooth_knot_distance",
                                    "smooth_knot_distance_power",
                                    "smooth_deg",
                                    "smooth_diff",
                                    "gamma");

  Rcpp::NumericMatrix chosen_parameters = Rcpp::wrap(chosen_params);
  Rcpp::colnames(chosen_parameters) = param_names;

  object["chosen_parameters"] = chosen_parameters;
  object["opt_index"] = opt_index + 1;
  object["predictions"] = predictions;
  object["weights"] = weights;
  object["past_performance"] = past_performance;
  object["forecaster_loss"] = loss_for;
  object["experts_loss"] = loss_exp;

  return object;
}