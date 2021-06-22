// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>

#include <loss.h>
#include <misc.h>
#include <splines.h>
#include <oracle.h>

using namespace arma;

//' @template function_batch
//'
//' @template param_y
//' @template param_experts
//' @template param_tau
//' @template param_affine
//' @template param_positive
//' @template param_intercept
//' @template param_debias
//' @param initial_window Defines the size of the initial estimaton window.
//' @param expanding_window Defines wether an expanding window or a rolling window shall be used for batch optimization. Defaults to TRUE.
//' @template param_loss_function
//' @template param_loss_parameter
//' @template param_lambda
//' @template param_forget
//' @template param_forget_performance
//' @template param_fixed_share
//' @param gamma to be removed
//' @template param_ndiff
//' @template param_deg
//' @template param_knot_distance
//' @template param_knot_distance_power
//' @template param_trace
//' @template param_lead_time
//' @template param_allow_quantile_crossing
//' @usage batch(y, experts, tau, affine = FALSE, positive = FALSE, intercept = FALSE,
//' debias = TRUE, initial_window = 30, expanding_window = TRUE,
//' loss_function = "quantile", loss_parameter = 1, lambda = -Inf,
//' forget = 0, forget_performance = 0, fixed_share = 0, gamma = 1, ndiff = 1, deg = 3,
//' knot_distance = 0.025, knot_distance_power = 1, trace = TRUE, lead_time = 0,
//' allow_quantile_crossing = FALSE)
//' @export
// [[Rcpp::export]]
Rcpp::List batch(
    mat &y,
    cube &experts,
    Rcpp::NumericVector tau = Rcpp::NumericVector::create(),
    const bool &affine = false,
    const bool &positive = false,
    const bool &intercept = false,
    const bool &debias = true,
    int initial_window = 30,
    const bool expanding_window = true,
    const std::string loss_function = "quantile",
    const double &loss_parameter = 1,
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(),
    Rcpp::NumericVector forget = Rcpp::NumericVector::create(),
    const double &forget_performance = 0,
    Rcpp::NumericVector fixed_share = Rcpp::NumericVector::create(),
    Rcpp::NumericVector gamma = Rcpp::NumericVector::create(),
    Rcpp::NumericVector ndiff = Rcpp::NumericVector::create(),
    Rcpp::NumericVector deg = Rcpp::NumericVector::create(),
    Rcpp::NumericVector knot_distance = Rcpp::NumericVector::create(),
    Rcpp::NumericVector knot_distance_power = Rcpp::NumericVector::create(),
    const bool trace = true,
    const int &lead_time = 0,
    bool allow_quantile_crossing = false,
    Rcpp::NumericVector soft_threshold = Rcpp::NumericVector::create(),
    Rcpp::NumericVector hard_threshold = Rcpp::NumericVector::create())
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

    if (T <= initial_window)
        Rcpp::stop("Initial estimation window exceeds input data.");

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

    // Set default values
    vec lambda_vec = set_default(lambda, -datum::inf);
    vec forget_vec = set_default(forget, 0);
    vec fixed_share_vec = set_default(fixed_share, 0);
    vec gamma_vec = set_default(gamma, 1);
    vec knot_distance_vec = set_default(knot_distance, 0.025);
    vec deg_vec = set_default(deg, 3);
    vec diff_vec = set_default(ndiff, 1.5);
    vec knots_asym_vec = set_default(knot_distance_power, 1);
    vec threshold_soft_vec = set_default(soft_threshold, -datum::inf);
    vec threshold_hard_vec = set_default(hard_threshold, -datum::inf);

    // Init parametergrid
    mat param_grid = get_combinations(lambda_vec, forget_vec);     // Index 0 & 1
    param_grid = get_combinations(param_grid, fixed_share_vec);    // Index 2
    param_grid = get_combinations(param_grid, gamma_vec);          // Index 3
    param_grid = get_combinations(param_grid, knot_distance_vec);  // Index 4
    param_grid = get_combinations(param_grid, deg_vec);            // Index 5
    param_grid = get_combinations(param_grid, diff_vec);           // Index 6
    param_grid = get_combinations(param_grid, knots_asym_vec);     // Index 7
    param_grid = get_combinations(param_grid, threshold_soft_vec); // Index 8
    param_grid = get_combinations(param_grid, threshold_hard_vec); // Index 9

    const int X = param_grid.n_rows;
    mat chosen_params(T, param_grid.n_cols);
    vec opt_index(T + 1, fill::zeros);
    cube past_performance(T, P, X, fill::zeros);
    vec tmp_performance(X, fill::zeros);
    vec cum_performance(X, fill::zeros);
    Progress prog((T - initial_window) * X + X, trace);

    // Populate uniform weights

    // Init object holding temp. weights, resp. ex-ante
    cube w_post(P, K, X);
    w_post.fill(1 / double(K - intercept * debias));

    if (intercept && debias)
        w_post.col(0).fill(0);

    // Weights output
    cube weights(T + 1, P, K, fill::zeros);

    mat experts_mat(P, K);

    cube predictions_post(T, P, X);
    mat predictions_temp(P, K);
    mat predictions_final(T + T_E_Y, P, fill::zeros);

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

    // Used for rolling window
    int start = 0;

    // Predictions at t < lead_time + initial_window  using initial weights
    for (unsigned int t = 0; t < initial_window + lead_time; t++)
    {
        // Save final weights w_post
        weights.row(t) = w_post.slice(opt_index(t));

        // Store expert predictions temporarily
        experts_mat = experts.row(t);

        // Forecasters prediction (ex-post)
        predictions_temp = sum(w_post.each_slice() % experts_mat, 1);
        predictions_post.row(t) = predictions_temp;

        // Final prediction
        predictions_final.row(t) =
            vectorise(predictions_post(span(t), span::all, span(opt_index(t)))).t();

        past_performance.row(t).fill(datum::nan);
    }

    for (unsigned int t = (0 + initial_window + lead_time); t < T; t++)
    {
        if (!expanding_window)
        {
            start += 1;
        }

        // Save final weights w_post
        weights.row(t) = w_post.slice(opt_index(t));

        // Store expert predictions temporarily
        experts_mat = experts.row(t);

        // Forecasters prediction (ex-post)
        predictions_temp = sum(w_post.each_slice() % experts_mat, 1);
        predictions_post.row(t) = predictions_temp;

        // Final prediction
        predictions_final.row(t) =
            vectorise(predictions_post(span(t), span::all, span(opt_index(t)))).t();

        for (unsigned int x = 0; x < X; x++)
        {
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

                // optim_weights()
                mat experts_tmp = experts.tube(span(start, t - lead_time), span(p, p));

                w_post(span(p), span::all, span(x)) = optimize_weights(
                    y(span(start, t - lead_time), p),
                    experts_tmp,
                    affine,
                    positive,
                    intercept,
                    debias,
                    loss_function,
                    tau_vec(p),
                    param_grid(x, 1), // Forget
                    loss_parameter);

                R_CheckUserInterrupt();
            }

            // Apply thresholds

            for (double &e : w_post(span::all, span(intercept * debias, K - 1), span(x)))
            {
                threshold_soft(e, param_grid(x, 8));
            }

            for (double &e : w_post(span::all, span(intercept * debias, K - 1), span(x)))
            {
                threshold_hard(e, param_grid(x, 9));
            }

            // Smoothing
            if (param_grid(x, 0) != -datum::inf)
            {
                for (unsigned int k = 0; k < K; k++)
                {
                    w_post(span::all, span(k), span(x)) =
                        hat_mats(x) *
                        vectorise(w_post(span::all, span(k), span(x)));
                }
            }

            for (unsigned int p = 0; p < P; p++)
            {
                //Add fixed_share
                w_post(span(p), span(intercept * debias, K - 1), span(x)) =
                    (1 - param_grid(x, 2)) * w_post(span(p), span(intercept * debias, K - 1), span(x)) +
                    (param_grid(x, 2) / (K - intercept * debias));

                // Ensure constraints are met
                if (positive)
                {
                    w_post(span(p), span(intercept * debias, K - 1), span(x)) = pmax_arma(w_post(span(p), span(intercept * debias, K - 1), span(x)), 0);
                }

                if (affine)
                {
                    w_post(span(p), span(intercept * debias, K - 1), span(x)) /= accu(w_post(span(p), span(intercept * debias, K - 1), span(x)));
                }
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

    // Set unused values to NA
    chosen_params.rows(0, lead_time + initial_window - 1).fill(datum::nan);

    // 1-Indexing for R-Output
    opt_index = opt_index + 1;

    Rcpp::DataFrame opt_params_df = Rcpp::DataFrame::create(
        Rcpp::Named("lambda") = chosen_params.col(0),
        Rcpp::Named("forget") = chosen_params.col(1),
        Rcpp::Named("fixed_share") = chosen_params.col(2),
        Rcpp::Named("gamma") = chosen_params.col(3),
        Rcpp::Named("knot_distance") = chosen_params.col(4),
        Rcpp::Named("deg") = chosen_params.col(5),
        Rcpp::Named("diff") = chosen_params.col(6),
        Rcpp::Named("knot_distance_power") = chosen_params.col(7));

    Rcpp::DataFrame parametergrid = Rcpp::DataFrame::create(
        Rcpp::Named("lambda") = param_grid.col(0),
        Rcpp::Named("forget") = param_grid.col(1),
        Rcpp::Named("fixed_share") = param_grid.col(2),
        Rcpp::Named("gamma") = param_grid.col(3),
        Rcpp::Named("knot_distance") = param_grid.col(4),
        Rcpp::Named("deg") = param_grid.col(5),
        Rcpp::Named("diff") = param_grid.col(6),
        Rcpp::Named("knot_distance_power") = param_grid.col(7));

    return Rcpp::List::create(
        Rcpp::Named("predictions") = predictions_final,
        Rcpp::Named("weights") = weights,
        Rcpp::Named("forecaster_loss") = loss_for,
        Rcpp::Named("experts_loss") = loss_exp,
        Rcpp::Named("past_perf_wrt_params") = past_performance,
        Rcpp::Named("chosen_parameters") = opt_params_df,
        Rcpp::Named("parametergrid") = parametergrid,
        Rcpp::Named("opt_index") = opt_index);
}