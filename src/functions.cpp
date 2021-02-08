// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <string>
#include <progress.hpp>
using namespace arma;

mat bbase(const vec &x, const int &nseg, const int &deg)
{

  double xl = min(x);
  double xr = max(x);

  double dx = (xr - xl) / nseg;
  int K = round(((xr + deg * dx) - (xl - deg * dx)) / dx + 1);
  vec knots = linspace(xl - deg * dx, xr + deg * dx, K);
  mat P(x.n_rows, knots.n_rows);

  for (unsigned int row_ = 0; row_ < P.n_rows; row_++)
  {
    for (unsigned int col_ = 0; col_ < P.n_cols; col_++)
    {
      P(row_, col_) = pow((x(row_) - knots(col_)), deg) * (x(row_) > knots(col_));
    }
  }

  mat D(knots.n_rows, knots.n_rows, fill::eye);
  D = diff(D, deg + 1);
  D = D / (exp(lgamma(deg + 1)) * pow(dx, deg));
  mat B = pow(-1, deg + 1) * P * D.t();

  return B;
}

vec pmax_arma(const vec &x,
              const double &bound)
{
  size_t n = x.n_elem;
  vec out = x;
  for (size_t i = 0; i < n; i++)
  {
    if (out(i) < bound)
      out(i) = bound;
  }
  return out;
}

vec pmin_arma(const vec &x,
              const double &bound)
{
  size_t n = x.n_elem;
  vec out = x;
  for (size_t i = 0; i < n; i++)
  {
    if (out(i) > bound)
      out(i) = bound;
  }
  return out;
}

double loss(const double &y,
            const double &x,
            const double &pred = 0,
            const std::string method = "ql",
            const double &tau = 0.5,
            const double &alpha = 2,
            const bool &gradient = true)
{
  double loss;

  if (gradient)
  {
    loss = ((y < pred) - tau) * x;
  }
  else
  {
    loss = ((y < x) - tau) * (x - y);
  }
  return loss;
}

mat get_combinations(const mat &x, const vec &y)
{

  int i = 0;

  mat grid(x.n_rows * y.size(), x.n_cols + 1);

  for (unsigned int y_ = 0; y_ < y.size(); y_++)
  {
    for (unsigned int x_ = 0; x_ < x.n_rows; x_++)
    {
      grid(i, span(0, x.n_cols - 1)) = x.row(x_);
      grid(i, x.n_cols) = y(y_);
      i += 1;
    }
  }
  return (grid);
}

vec set_default(const vec &input,
                const double &value)
{
  vec output = input;
  if (output.size() == 0)
  {
    output.set_size(1);
    output(0) = value;
  }

  return output;
}

//' Probabilistic Forecast Combination - ProFoC
//'
//' Returns predictions and weights calculated by online-learning algorithms
//' using CRPS Learning. By default, the weights are calculated by
//' gradient based bernstein online aggregation (BOAG).
//'
//' @param y A numeric vector of realizations.
//' @param experts A an array of predictions with dimension
//' (Observations, Quantiles, Experts).
//' @param tau A numeric vector of probabilities.
//' corresponding to the columns of experts.
//' @param ex_post_smooth Determines if smoothing is during or after
//' online-learning. If true, contemporary weights are not affected
//' but output weights are. If false (default) smoothed weights are
//' also by the algorithm.
//' @param ex_post_fs Analogous to ex_post_smooth: shall a fixed-share
//' be added during (FALSE) or after online-learning (TRUE).
//' @param lambda Penalization parameter used in the smoothing Step.
//' -Inf causes the smoothing step to be skipped (default).
//' @param method One of "boa", "ml_poly" or "ewa".
//' @param method_var Allows to calculate slight variations of the BOA
//' algorithm
//' @param forget Share of past Regret not to be considered, resp. to be
//' forgotten in every iteration of the algorithm. Defaults to 0.
//' @param fixed_share Amount of fixed share to be added to the weights.
//' Defaults to 0. 1 leads to uniform weights.
//' @param gamma Scaling parameter for the learning rate.
//' @param ndiff Degree of the differencing operator in the smoothing equation. 1 (default) leads to shrikage towards a constant.
//' @param deg Degree of the B-Spine basis functions.
//' @param rel_nseg determines how many basis functions are created in the smoothing step. This parameter defaults to 0.5 leading to 0.5*length(tau) to be used as number of knots. 1 Leads to as many knots as tau is long. 0 leads to one single knot.
//' @param gradient Determines if a linearized version of the loss is used.
//' @param loss_array User specified loss array. If specified, the loss will not be calculated by profoc.
//' @param regret_array User specified regret array. If specifiec, the regret will not be calculated by profoc.
//' @usage profoc(y, experts, tau, ex_post_smooth = FALSE, ex_post_fs = FALSE,
//' lambda = -Inf, method = "boa", method_var = "A", forget = 0,
//' fixed_share = 0, gamma = 1, ndiff = 1, deg = 3, rel_nseg = 0.5,
//' gradient = TRUE, loss_array = NULL, regret_array = NULL)
//' @return Profoc can tune various parameters automatically based on
//' the past loss. For this, lambda, forget, fixed_share, gamma, ndiff,
//' deg and rel_nseg can be specified as numeric vectors containing
//' parameters to consider. Profoc will automatically try all possible
//' combinations of values provide.
//' @export
// [[Rcpp::export]]
Rcpp::List profoc(
    mat &y,
    const cube &experts,
    Rcpp::NumericVector tau = Rcpp::NumericVector::create(),
    const bool &ex_post_smooth = false,
    const bool &ex_post_fs = false,
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(),
    const std::string method = "boa",
    const std::string method_var = "A",
    Rcpp::NumericVector forget = Rcpp::NumericVector::create(), Rcpp::NumericVector fixed_share = Rcpp::NumericVector::create(), Rcpp::NumericVector gamma = Rcpp::NumericVector::create(), Rcpp::NumericVector ndiff = Rcpp::NumericVector::create(), Rcpp::NumericVector deg = Rcpp::NumericVector::create(), Rcpp::NumericVector rel_nseg = Rcpp::NumericVector::create(),
    const bool &gradient = true,
    Rcpp::NumericVector loss_array = Rcpp::NumericVector::create(), Rcpp::NumericVector regret_array = Rcpp::NumericVector::create())
{

  // Indexing Convention -> (T, P, K, X)
  // T number of observations
  // P lengths of probability Grid
  // K number of experts
  // X number of parameter combinations to consider

  // Object Dimensions
  const int T = experts.n_rows;
  const int P = experts.n_cols;
  const int K = experts.n_slices;

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

  // Set default values
  vec lambda_vec = set_default(lambda, -datum::inf);
  vec forget_vec = set_default(forget, 0);
  vec fixed_share_vec = set_default(fixed_share, 0);
  vec gamma_vec = set_default(gamma, 1);
  vec rel_nseg_vec = set_default(rel_nseg, 0.5);
  vec deg_vec = set_default(deg, 3);
  vec diff_vec = set_default(ndiff, 1);

  // Init parametergrid
  mat param_grid = get_combinations(lambda_vec, forget_vec);
  param_grid = get_combinations(param_grid, fixed_share_vec);
  param_grid = get_combinations(param_grid, gamma_vec);
  param_grid = get_combinations(param_grid, rel_nseg_vec);
  param_grid = get_combinations(param_grid, deg_vec);
  param_grid = get_combinations(param_grid, diff_vec);

  const int X = param_grid.n_rows;
  mat chosen_params(T, param_grid.n_cols);
  int opt_index = 0;
  vec opt_index_(T);
  cube past_performance(T, P, X, fill::zeros);

  Progress prog(T * X + X, true);

  cube weights(T + 1, P, K, fill::zeros);
  cube w(P, K, X);
  w = w.fill(1 / double(K));
  cube w_temp(P, K, X, fill::ones);
  w_temp = w_temp.fill(1 / double(K));
  vec w0(K);
  w0 = w0.fill(1 / double(K));
  cube R(P, K, X, fill::zeros);
  cube R_reg(R);
  cube eta(P, K, X, fill::zeros);
  if (method == "ml_poly")
    eta.fill(exp(350));
  cube V(P, K, X, fill::zeros);
  cube E(P, K, X, fill::zeros);
  cube k(P, K, X, fill::zeros);
  k = k.fill(-699);
  vec experts_vec(K);
  mat experts_mat(P, K);
  vec forecasters_pred(1);
  mat predictions(T, P);
  vec lpred(1);
  vec lexp(K);
  vec r(K);
  vec r_reg(K);
  mat mean_pb_loss(param_grid.n_rows, P);
  vec crps_approx(param_grid.n_rows);
  double y_tilde;
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
  vec spline_basis_x = regspace(1, P) / (P + 1);

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
          param_grid(x, 6) == param_grid(x - 1, 6))
      {
        hat_mats(x) = hat_mats(x - 1);
      }
      else
      {
        // Number of segments
        int nseg = std::min(std::max(double(1), ceil(P * param_grid(x, 4))), double(P));

        // Create bspline basis
        mat B = bbase(spline_basis_x, nseg, param_grid(x, 5));

        // Create roughness measure
        mat D(B.n_cols, B.n_cols, fill::eye);
        D = diff(D, param_grid(x, 6));

        hat_mats(x) = B * inv(B.t() * B + param_grid(x, 0) * D.t() * D) * B.t();
      }
      R_CheckUserInterrupt();
      prog.increment(); // Update progress
    }
  }

  for (unsigned int t = 0; t < T; t++)
  {

    // Save Weights and Prediction
    // Note that w_temp != w in ex post setting and w_temp = w otherwise
    experts_mat = experts.row(t);
    weights.row(t) = w_temp.slice(opt_index);
    predictions.row(t) = sum(w_temp.slice(opt_index) % experts_mat, 1).t();

    for (unsigned int x = 0; x < X; x++)
    {

      for (unsigned int p = 0; p < P; p++)
      {

        // Forecasters prediction w.r.t. specific parameters (ex-post)
        y_tilde = sum(
            vectorise(w_temp(span(p), span::all, span(x))) %
            vectorise(experts(span(t), span(p), span::all)));
        experts_vec = experts(span(t), span(p), span::all);
        // Forecasters prediction w.r.t. specific parameters (ex-ante)
        forecasters_pred = experts_vec.t() *
                           vectorise(w(span(p), span::all, span(x)));

        for (unsigned int k = 0; k < K; k++)
        {
          past_performance(t, p, x) = loss(y(t, p),
                                           y_tilde,
                                           9999,       // where to evaluate gradient
                                           "ql",       // method
                                           tau_vec(p), // tau_vec
                                           2,          // alpha
                                           false);
          lexp(k) = loss(y(t, p),
                         experts_vec(k),
                         forecasters_pred(0), // where to evaluate gradient
                         "ql",                // method
                         tau_vec(p),          // tau_vec
                         2,                   // alpha
                         gradient);
        }

        lpred = loss(y(t, p),
                     forecasters_pred(0),
                     forecasters_pred(0), // where to evaluate gradient
                     "ql",                // method
                     tau_vec(p),          // tau_vec
                     2,                   // alpha
                     gradient);

        if (loss_array.size() != 0)
        {
          lexp = loss_cube(span(t), span(p), span::all);
        }

        if (regret_array.size() == 0)
        {
          r = lpred - lexp.each_row();
        }
        else
        {
          r = regret_cube(span(t), span(p), span::all);
        }

        if (method == "ewa")
        {
          // Update the cumulative regret used by eta
          R(span(p), span::all, span(x)) = vectorise(R(span(p), span::all, span(x))) * (1 - param_grid(x, 1)) + r;
          w(span(p), span::all, span(x)) = exp(param_grid(x, 3) * vectorise(R(span(p), span::all, span(x))));
          w(span(p), span::all, span(x)) = w(span(p), span::all, span(x)) / accu(w(span(p), span::all, span(x)));
        }
        else if (method == "ml_poly")
        {
          // Update the cumulative regret used by ML_Poly
          R(span(p), span::all, span(x)) = vectorise(R(span(p), span::all, span(x))) * (1 - param_grid(x, 1)) + r;

          // Update the learning rate
          eta(span(p), span::all, span(x)) = 1 / (1 / vectorise(eta(span(p), span::all, span(x))) + square(r));

          // param_grid(x, 3) = gamma
          w(span(p), span::all, span(x)) = param_grid(x, 3) * vectorise(eta(span(p), span::all, span(x))) % pmax_arma(vectorise(R(span(p), span::all, span(x))), exp(-700));
          w(span(p), span::all, span(x)) = w(span(p), span::all, span(x)) / accu(w(span(p), span::all, span(x)));
        }
        else if (method == "boa")
        {
          // Opposite of instantaneous regret
          r = -r;

          V(span(p), span::all, span(x)) = vectorise(V(span(p), span::all, span(x))) * (1 - param_grid(x, 1)) + square(r);

          // Update the bounds and learning rates
          E(span(p), span::all, span(x)) = max(vectorise(E(span(p), span::all, span(x))) * (1 - param_grid(x, 1)), abs(r));

          eta(span(p), span::all, span(x)) =
              pmin_arma(
                  min(1 / (2 * vectorise(E(span(p), span::all, span(x)))),
                      sqrt(log(K) / vectorise(V(span(p), span::all, span(x))))),
                  exp(350));

          // Instantaneous regularized regret
          r_reg = r + vectorise(eta(span(p), span::all, span(x))) % square(r);

          // Update cumulative regularized regret
          if (method_var.find('A') != std::string::npos)
          {
            R_reg(span(p), span::all, span(x)) = vectorise(R_reg(span(p), span::all, span(x))) * (1 - param_grid(x, 1)) + 0.5 * (r_reg + conv_to<colvec>::from(vectorise(eta(span(p), span::all, span(x))) % r > 0.5) % (2 * vectorise(E(span(p), span::all, span(x)))));
          }
          else
          {
            // Gaillard
            R_reg(span(p), span::all, span(x)) = vectorise(R_reg(span(p), span::all, span(x))) * (1 - param_grid(x, 1)) + r_reg;
          }

          // Update weights and truncate if necessary
          if (method_var.find('C') != std::string::npos)
          {
            // Wintenberger
            w(span(p), span::all, span(x)) = param_grid(x, 3) * vectorise(eta(span(p), span::all, span(x))) % exp(-param_grid(x, 3) * vectorise(eta(span(p), span::all, span(x))) % vectorise(R_reg(span(p), span::all, span(x)))) % w0;
            w(span(p), span::all, span(x)) = w(span(p), span::all, span(x)) / mean(param_grid(x, 3) * vectorise(eta(span(p), span::all, span(x))) % exp(-param_grid(x, 3) * vectorise(eta(span(p), span::all, span(x))) % vectorise(R_reg(span(p), span::all, span(x)))));
          }
          else
          {
            // Gaillard
            w(span(p), span::all, span(x)) = exp(-param_grid(x, 3) * vectorise(eta(span(p), span::all, span(x))) % vectorise(R_reg(span(p), span::all, span(x))));
            w(span(p), span::all, span(x)) = pmin_arma(pmax_arma(vectorise(w(span(p), span::all, span(x))), exp(-700)), exp(700));
            w(span(p), span::all, span(x)) = w(span(p), span::all, span(x)) / accu(w(span(p), span::all, span(x)));
          }
        }
        else
        {
          Rcpp::stop("Choose 'boa', 'ml_poly' or 'ewa' as method.");
        }
      }

      w_temp.slice(x) = w.slice(x);

      // Smoothing Step
      if (param_grid(x, 0) != -datum::inf)
      {
        for (unsigned int k = 0; k < K; k++)
        {
          if (!ex_post_smooth)
          {
            w(span::all, span(k), span(x)) = hat_mats(x) * vectorise(w(span::all, span(k), span(x)));
            // Truncate to zero
            w(span::all, span(k), span(x)) = pmax_arma(w(span::all, span(k), span(x)), exp(-700));
            w_temp(span::all, span(k), span(x)) = w(span::all, span(k), span(x));
          }
          else
          {
            w_temp(span::all, span(k), span(x)) = hat_mats(x) * vectorise(w(span::all, span(k), span(x)));
            w_temp(span::all, span(k), span(x)) = pmax_arma(w_temp(span::all, span(k), span(x)), exp(-700));
          }
        }
      }

      //Add fixed_share
      for (unsigned int p = 0; p < P; p++)
      {
        if (ex_post_fs) // Fixed share only added to w_temp
        {
          w_temp(span(p), span::all, span(x)) = (1 - param_grid(x, 2)) * (w_temp(span(p), span::all, span(x)) / accu(w_temp(span(p), span::all, span(x)))) + (param_grid(x, 2) / K);
        }
        else if (ex_post_smooth && !ex_post_fs)
        {
          // Add fixed_share to wtemp
          w_temp(span(p), span::all, span(x)) = (1 - param_grid(x, 2)) * (w_temp(span(p), span::all, span(x)) / accu(w_temp(span(p), span::all, span(x)))) + (param_grid(x, 2) / K);
          // Add fixed_share to w
          w(span(p), span::all, span(x)) = (1 - param_grid(x, 2)) * (w(span(p), span::all, span(x)) / accu(w(span(p), span::all, span(x)))) + (param_grid(x, 2) / K);
        }
        else if (!ex_post_smooth && !ex_post_fs)
        {
          w(span(p), span::all, span(x)) = (1 - param_grid(x, 2)) * (w(span(p), span::all, span(x)) / accu(w(span(p), span::all, span(x)))) + (param_grid(x, 2) / K);
          w_temp(span(p), span::all, span(x)) = w(span(p), span::all, span(x));
        }
      }
      prog.increment(); // Update progress
      R_CheckUserInterrupt();
    }

    // Sum past_performance in each slice
    crps_approx = sum(sum(past_performance, 0), 1);
    opt_index = crps_approx.index_min();
    opt_index_(t) = opt_index + 1;
    chosen_params.row(t) = param_grid.row(opt_index);
  }

  // Save Weights and Prediction
  weights.row(T) = w_temp.slice(opt_index);

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
                 9999,       // where to evaluate the gradient
                 "ql",       // method
                 tau_vec(p), // tau_vec
                 2,          // alpha
                 false);     // gradient
      }
      loss_for(t, p) = loss(y(t, p),
                            predictions(t, p),
                            9999,       // where to evaluate the gradient
                            "ql",       // method
                            tau_vec(p), // tau_vec
                            2,          // alpha
                            false);     // gradient;
    }
  }

  Rcpp::DataFrame opt_params_df = Rcpp::DataFrame::create(
      Rcpp::Named("lambda") = chosen_params.col(0),
      Rcpp::Named("forget") = chosen_params.col(1),
      Rcpp::Named("fixed_share") = chosen_params.col(2),
      Rcpp::Named("gamma") = chosen_params.col(3),
      Rcpp::Named("rel_nseg") = chosen_params.col(4),
      Rcpp::Named("deg") = chosen_params.col(5),
      Rcpp::Named("diff") = chosen_params.col(6));
  ;

  Rcpp::DataFrame parametergrid = Rcpp::DataFrame::create(
      Rcpp::Named("lambda") = param_grid.col(0),
      Rcpp::Named("forget") = param_grid.col(1),
      Rcpp::Named("fixed_share") = param_grid.col(2),
      Rcpp::Named("gamma") = param_grid.col(3),
      Rcpp::Named("rel_nseg") = param_grid.col(4),
      Rcpp::Named("deg") = param_grid.col(5),
      Rcpp::Named("diff") = param_grid.col(6));
  ;

  return Rcpp::List::create(
      Rcpp::Named("predictions") = predictions,
      Rcpp::Named("weights") = weights,
      Rcpp::Named("pb_loss_forecaster") = loss_for,
      Rcpp::Named("pb_loss_experts") = loss_exp,
      Rcpp::Named("past_perf_wrt_params") = past_performance,
      Rcpp::Named("chosen_parameters") = opt_params_df,
      Rcpp::Named("parametergrid") = parametergrid,
      Rcpp::Named("opt_index") = opt_index_);
}

//' Spline Fit - Fit Penalized B-Splines
//'
//' Returns fitted values.
//'
//' @param y The respones variable. Must be a numeric vector.
//' @param x Explanatory variable. Must be a numeric vector.
//' @param lambda The penalization parameter.
//' @param ndiff Degree of the difference operator to use when calculating
//' the B-Spline basis. A value of 1 corresponds to shrinkage towards
//' a constant, 2 corresponds to shrinkage towards a linear relationship
//' between x and y.
//' @param deg the degree of the basis functions.
//' @param rel_nseg The number of basis functions are calculated relative to the
//' length of y. A value of 0 corresponds to only 1 knot. A value of 1
//' corresponds to length(y) knots.
//' @export
// [[Rcpp::export]]
vec spline_fit(const vec &y, const vec &x, const double &lambda = 1, const int &ndiff = 1, const int &deg = 3, const double &rel_nseg = 0.1)
{

  // Number of segments
  int nseg = std::max(double(1), ceil(y.n_rows * rel_nseg));

  // Create bspline basis
  mat B = bbase(x, nseg, deg);

  // Create roughness measure
  mat D(B.n_cols, B.n_cols, fill::eye);
  D = diff(D, ndiff);

  // Calculate fitted values
  vec y_hat = B * inv(B.t() * B + lambda * D.t() * D) * B.t() * y;

  return y_hat;
}