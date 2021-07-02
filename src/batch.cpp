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
//' @template param_lead_time
//'
//' @param initial_window Defines the size of the initial estimaton window.
//' @param expanding_window Defines wether an expanding window or a rolling window shall be used for batch optimization. Defaults to TRUE.
//'
//' @template param_loss_function
//' @template param_loss_parameter
//'
//' @template param_basis_knot_distance_batch
//' @template param_basis_knot_distance_power
//' @template param_basis_deg_batch
//'
//' @template param_forget
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
//' @template param_parametergrid_max_combinations
//' @template param_forget_past_performance
//'
//' @template param_allow_quantile_crossing
//'
//' @template param_trace
//'
//' @usage batch(y,
//' experts,
//' tau,
//' affine = FALSE,
//' positive = FALSE,
//' intercept = FALSE,
//' debias = TRUE,
//' lead_time = 0,
//' initial_window = 30,
//' expanding_window = TRUE,
//' loss_function = "quantile",
//' loss_parameter = 1,
//' basis_knot_distance = 0.01,
//' basis_knot_distance_power = 1,
//' basis_deg = 1,
//' forget = 0,
//' soft_threshold = -Inf,
//' hard_threshold = -Inf,
//' fixed_share = 0,
//' smooth_lambda = -Inf,
//' smooth_knot_distance = 0.01,
//' smooth_knot_distance_power = 1,
//' smooth_deg = 1,
//' smooth_ndiff = 1,
//' parametergrid_max_combinations = 100,
//' forget_past_performance = 0,
//' allow_quantile_crossing = FALSE,
//' trace = TRUE
//' )
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
    const int &lead_time = 0,
    int initial_window = 30,
    const bool expanding_window = true,
    const std::string loss_function = "quantile",
    const double &loss_parameter = 1,
    Rcpp::NumericVector basis_knot_distance = Rcpp::NumericVector::create(),
    Rcpp::NumericVector basis_knot_distance_power = Rcpp::NumericVector::create(1),
    Rcpp::NumericVector basis_deg = Rcpp::NumericVector::create(1),
    Rcpp::NumericVector forget = Rcpp::NumericVector::create(0),
    Rcpp::NumericVector soft_threshold = Rcpp::NumericVector::create(-1 / 0),
    Rcpp::NumericVector hard_threshold = Rcpp::NumericVector::create(-1 / 0),
    Rcpp::NumericVector fixed_share = Rcpp::NumericVector::create(0),
    Rcpp::NumericVector smooth_lambda = Rcpp::NumericVector::create(-1 / 0),
    Rcpp::NumericVector smooth_knot_distance = Rcpp::NumericVector::create(),
    Rcpp::NumericVector smooth_knot_distance_power = Rcpp::NumericVector::create(),
    Rcpp::NumericVector smooth_deg = Rcpp::NumericVector::create(),
    Rcpp::NumericVector smooth_ndiff = Rcpp::NumericVector::create(1.5),
    const int parametergrid_max_combinations = 100,
    const double &forget_past_performance = 0,
    bool allow_quantile_crossing = false,
    const bool trace = true)
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

    vec basis_knot_distance_vec = set_default(basis_knot_distance, 1 / (double(P) + 1));

    bool deg_inheritance = false;
    if (smooth_deg.size() == 0)
        deg_inheritance = true;

    bool knot_distance_inheritance = false;
    if (smooth_knot_distance.size() == 0)
        knot_distance_inheritance = true;

    bool knot_distance_power_inheritance = false;
    if (smooth_knot_distance_power.size() == 0)
        knot_distance_power_inheritance = true;

    // Init parametergrid
    mat param_grid =
        get_combinations(basis_knot_distance_vec,                                                              // Index 0
                         basis_knot_distance_power);                                                           // Index 1
    param_grid = get_combinations(param_grid, basis_deg);                                                      // index 2
    param_grid = get_combinations(param_grid, forget);                                                         // index 3
    param_grid = get_combinations(param_grid, soft_threshold);                                                 // index 4
    param_grid = get_combinations(param_grid, hard_threshold);                                                 // index 5
    param_grid = get_combinations(param_grid, fixed_share);                                                    // index 6
    param_grid = get_combinations(param_grid, smooth_lambda);                                                  // index 7
    param_grid = get_combinations(param_grid, smooth_knot_distance, knot_distance_inheritance, 0);             // Index 8
    param_grid = get_combinations(param_grid, smooth_knot_distance_power, knot_distance_power_inheritance, 1); // Index 9
    param_grid = get_combinations(param_grid, smooth_deg, deg_inheritance, 2);                                 // Index 10
    param_grid = get_combinations(param_grid, smooth_ndiff);                                                   // Index 11

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

    field<sp_mat> basis_mats(X);
    field<mat> beta(X);

    // Init hat matrix field
    for (unsigned int x = 0; x < X; x++)
    {
        mat basis;

        // In second step: skip if recalc is not necessary:
        if (x > 0 &&
            param_grid(x, 2) == param_grid(x - 1, 2) &&
            param_grid(x, 0) == param_grid(x - 1, 0) &&
            param_grid(x, 1) == param_grid(x - 1, 1))
        {
            basis_mats(x) = basis_mats(x - 1);
            basis = basis_mats(x);
        }
        else
        {
            basis = make_basis_matrix(spline_basis_x,
                                      param_grid(x, 0),  // kstep
                                      param_grid(x, 2),  // degree
                                      param_grid(x, 1)); // uneven grid

            basis_mats(x) = sp_mat(basis);
        }

        beta(x) = (w_post.slice(x).t() * pinv(basis).t()).t();

        R_CheckUserInterrupt();
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

        cube experts_tmp_cube = experts.rows(start, t - lead_time);

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

                if (basis_mats(x).is_diagmat())
                {
                    // optim_weights()
                    mat experts_tmp = experts_tmp_cube.col(p);

                    w_post(span(p), span::all, span(x)) = optimize_weights(
                        y(span(start, t - lead_time), p),
                        experts_tmp,
                        affine,
                        positive,
                        intercept,
                        debias,
                        loss_function,
                        tau_vec(p),
                        param_grid(x, 3), // Forget
                        loss_parameter);

                    R_CheckUserInterrupt();
                }
            }

            if (!basis_mats(x).is_diagmat())
            {
                beta(x) = optimize_betas(
                    y.rows(start, t - lead_time),
                    experts_tmp_cube,
                    affine,
                    positive,
                    intercept,
                    debias,
                    loss_function,
                    tau_vec,
                    param_grid(x, 3), // Forget
                    loss_parameter,
                    basis_mats(x),
                    beta(x));

                w_post.slice(x) = basis_mats(x) * beta(x);

                R_CheckUserInterrupt();
            }

            // Apply thresholds

            for (double &e : w_post(span::all, span(intercept * debias, K - 1), span(x)))
            {
                threshold_soft(e, param_grid(x, 4));
            }

            for (double &e : w_post(span::all, span(intercept * debias, K - 1), span(x)))
            {
                threshold_hard(e, param_grid(x, 5));
            }

            // Smoothing
            if (param_grid(x, 7) != -datum::inf)
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
                    (1 - param_grid(x, 6)) * w_post(span(p), span(intercept * debias, K - 1), span(x)) +
                    (param_grid(x, 6) / (K - intercept * debias));

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
        cum_performance = (1 - forget_past_performance) * cum_performance + tmp_performance;
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

    Rcpp::CharacterVector param_names =
        Rcpp::CharacterVector::create("basis_knot_distance",
                                      "basis_knot_distance_power",
                                      "basis_deg",
                                      "forget",
                                      "threshold_soft",
                                      "threshold_hard",
                                      "fixed_share",
                                      "smooth_lambda",
                                      "smooth_knot_distance",
                                      "smooth_knot_distance_power",
                                      "smooth_deg",
                                      "smooth_diff");

    Rcpp::NumericMatrix parametergrid = Rcpp::wrap(param_grid);
    Rcpp::colnames(parametergrid) = param_names;

    Rcpp::NumericMatrix chosen_parameters = Rcpp::wrap(chosen_params);
    Rcpp::colnames(chosen_parameters) = param_names;

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("predictions") = predictions_final,
        Rcpp::Named("weights") = weights,
        Rcpp::Named("forecaster_loss") = loss_for,
        Rcpp::Named("experts_loss") = loss_exp,
        Rcpp::Named("past_performance") = past_performance,
        Rcpp::Named("chosen_parameters") = chosen_parameters,
        Rcpp::Named("parametergrid") = parametergrid,
        Rcpp::Named("opt_index") = opt_index,
        Rcpp::Named("basis_matrices") = basis_mats);

    out.attr("class") = "profoc_batch";

    return out;
}