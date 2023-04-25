#include <RcppArmadillo.h>
#include <progress.hpp>
#include <thread>

#include "misc.h"
#include "loss.h"
#include "clock.h"
#include "conline.h"
#include "profoc_types.h"

using namespace arma;

// conline class was exposed via "profoc_types.h"
// So we can use it here as input and output if necessary

// Example:
// [[Rcpp::export]]
bool test_class_input(conline &obj)
{
    return obj.trace;
}

// [[Rcpp::export]]
conline test_class_output()
{
    conline instance;
    instance.trace = false;
    return instance;
}

// Constructor

// Getters

// Methods

void conline::set_defaults()
{
    // Expand tau if necessary
    if (tau.n_elem == 1)
    {
        tau.resize(P);
        tau.fill(tau(0));
    }

    // Initial params
    w0.ones(D, P, K);
    w0 /= K;

    R0.zeros(D, P, K);
}

void conline::set_grid_objects()
{
    opt_index.zeros(T + 1);
    if (save_past_performance)
    {
        past_performance.set_size(T);
    }
    else
    {
        past_performance.set_size(1);
    }
    tmp_performance.zeros(D, P);
    cum_performance.zeros(X);

    predictions.zeros(T + T_E_Y, D, P);
    if (save_predictions_grid)
    {
        predictions_grid.set_size(T + T_E_Y);
    }
    else
    {
        predictions_grid.set_size(1 + lead_time);
    }
    weights_tmp.set_size(X);
    weights.set_size(T + 1);

    loss_for.zeros(T, D, P);
    loss_exp.set_size(T);

    V.set_size(X);
    E.set_size(X);
    eta.set_size(X);
    R.set_size(X);
    beta.set_size(X);
    beta0field.set_size(X);

#define vr arma::vectorise

    for (unsigned int x = 0; x < X; x++)
    {
        clock.tick("init");
        unsigned int Pr = basis_pr(params["basis_pr_idx"](x) - 1).n_cols;
        unsigned int Dr = basis_mv(params["basis_mv_idx"](x) - 1).n_cols;

        // Learning parameters
        V(x).zeros(Dr, Pr, K);
        E(x).zeros(Dr, Pr, K);

        arma::cube eta_(Dr, Pr, K, fill::zeros);
        eta(x) = eta_;
        if (method == "ml_poly" || method == "potential")
        {
            eta_.fill(exp(350));
            eta(x) = eta_;
        }

        R(x).set_size(Dr, Pr, K);
        beta(x).set_size(Dr, Pr, K);
        weights_tmp(x).set_size(D, P, K);

        for (unsigned int d = 0; d < D; d++)
        {
            weights_tmp(x).row(d) = w0.row(d);
        }

        for (unsigned int k = 0; k < K; k++)
        {
            R(x).slice(k) = basis_mv(params["basis_mv_idx"](x) - 1).t() *
                            R0.slice(k) *
                            basis_pr(params["basis_pr_idx"](x) - 1);
            R(x).slice(k) = basis_mv(params["basis_mv_idx"](x) - 1).t() *
                            R0.slice(k) *
                            basis_pr(params["basis_pr_idx"](x) - 1);
            beta(x).slice(k) = pinv(
                                   mat(basis_mv(params["basis_mv_idx"](x) - 1))) *
                               w0.slice(k) *
                               pinv(mat(basis_pr(params["basis_pr_idx"](x) - 1))).t();
        }
        beta0field(x) = beta(x);
    }

    // Predictions at t < lead_time using initial weights
    for (unsigned int t = 0; t < lead_time; t++)
    {

        weights(t).set_size(D, P, K);

        if (save_past_performance)
        {
            past_performance(t).set_size(D, P, X);
            past_performance(t).fill(datum::nan);
        }
        predictions_grid(t).set_size(D, P, X);

        // Store predictions w.r.t. grid for time t
        cube tmp_preds_cube(D, P, X);

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
                vec tmp_preds_vec = sum(weights_temp % experts_mat, 1);

                tmp_preds_cube(span(d), span::all, span(x)) = tmp_preds_vec;
            }
        }
        predictions_grid(t) = tmp_preds_cube;
        // Final Prediction
        predictions.row(t) = tmp_preds_cube.slice(opt_index(t));
    }

    start = lead_time;
    clock.tock("init");
}

void conline::learn()
{
    Progress prog(T, trace);
    clock.tick("core");

    for (unsigned int tp = 0; tp < predictions_grid.n_rows; tp++)
    {
        predictions_grid(tp).set_size(D, P, X);
    }

    for (unsigned int t = start; t < T; t++)
    {

        weights(t).set_size(D, P, K);

        if (save_past_performance)
            past_performance(t).set_size(D, P, X);

        clock.tick("loss");

        // Store predictions w.r.t. grid for time t
        cube tmp_preds_cube(D, P, X);

        for (unsigned int d = 0; d < D; d++)
        {
            // Save final weights weights_tmp
            weights(t).row(d) = weights_tmp(opt_index(t)).row(d);

            // Store expert predictions temporarily
            mat experts_mat = experts(t).row(d);

            // Predictions using different parameter values
            for (unsigned int x = 0; x < X; x++)
            {

                mat weights_temp = weights_tmp(x).row(d);
                vec tmp_preds_vec = sum(weights_temp % experts_mat, 1);

                // Sort predictions if quantile_crossing is prohibited
                if (!allow_quantile_crossing)
                {
                    tmp_preds_vec = arma::sort(tmp_preds_vec, "ascend", 0);
                }
                tmp_preds_cube(span(d), span::all, span(x)) = tmp_preds_vec;
            }
        }

        // Define which observation to base the weight update on. That is the
        // most recent prediction unless lead_time is greater than 0.
        arma::uword predictions_tmp_idx;
        if (save_predictions_grid)
        {
            predictions_tmp_idx = t - lead_time;
            predictions_grid(t) = tmp_preds_cube;
        }
        else
        {
            // If save_predictions_grid is false, the first observation is
            // allways the one we base the weight update on. This is because
            // the size of predictions_grid is only 1 + lead_time.
            predictions_tmp_idx = 0;
            for (unsigned int tp = 0; tp < predictions_grid.n_rows - 1; tp++)
            {
                predictions_grid(tp) = predictions_grid(tp + 1);
            }
            predictions_grid(predictions_grid.n_rows - 1) = tmp_preds_cube;
        }

        // Final prediction
        predictions.row(t) = tmp_preds_cube.slice(opt_index(t));

        clock.tock("loss");
        for (unsigned int x = 0; x < X; x++)
        {
            clock.tick("regret");
            mat lexp_int(P, K); // Experts loss
            mat lexp_ext(P, K); // Experts loss
            mat lexp(P, K);     // Experts loss
            vec lfor(P);        // Forecasters loss
            cube regret_tmp(D, P, K);
            cube regret(basis_mv(params["basis_mv_idx"](x) - 1).n_cols,
                        basis_pr(params["basis_pr_idx"](x) - 1).n_cols,
                        K); // Dr x Pr x K

            for (unsigned int d = 0; d < D; d++)
            {
#pragma omp parallel for
                for (unsigned int p = 0; p < P; p++)
                {
                    tmp_performance(d, p) = loss(y(t, d),
                                                 predictions_grid(predictions_tmp_idx)(d, p, x),
                                                 9999,           // where evaluate loss_gradient
                                                 loss_function,  // method
                                                 tau(p),         // tau
                                                 loss_parameter, // alpha
                                                 false);

                    if (params["loss_share"](x) != 1)
                    {
                        for (unsigned int k = 0; k < K; k++)
                        {
                            lexp_int(p, k) = loss(y(t, d),
                                                  experts(t)(d, p, k),
                                                  predictions_grid(predictions_tmp_idx)(d, p, x), // where evaluate loss_gradient
                                                  loss_function,                                  // method
                                                  tau(p),                                         // tau
                                                  loss_parameter,                                 // alpha
                                                  loss_gradient);
                        }

                        if (params["loss_share"](x) == 0)
                        {
                            lexp.row(p) = lexp_int.row(p);
                        }
                        else
                        {
                            lexp_ext.row(p) = vr(loss_array(t).tube(d, p)).t();
                            lexp.row(p) = (1 - params["loss_share"](x)) * lexp_int.row(p) + params["loss_share"](x) * lexp_ext.row(p);
                        }
                    }
                    else
                    {
                        lexp_ext.row(p) = vr(loss_array(t).tube(d, p)).t();
                        lexp.row(p) = lexp_ext.row(p);
                    }
                    lfor(p) = loss(y(t, d),
                                   predictions_grid(predictions_tmp_idx)(d, p, x),
                                   predictions_grid(predictions_tmp_idx)(d, p, x), // where to evaluate loss_gradient
                                   loss_function,                                  // method
                                   tau(p),                                         // tau
                                   loss_parameter,                                 // alpha
                                   loss_gradient);
                }

                mat regret_int(P, K);
                mat regret_ext(P, K);

                if (params["regret_share"](x) != 1)
                {
                    regret_int = (lfor - lexp_int.each_col()).t();

                    if (params["regret_share"](x) == 0)
                    {
                        regret_tmp.row(d) = regret_int.t();
                    }
                    else
                    {
                        regret_ext = regret_array(t).row(d);
                        regret_ext = regret_ext.t();
                        regret_tmp.row(d) = ((1 - params["regret_share"](x)) * regret_int + params["regret_share"](x) * regret_ext).t();
                    }
                }
                else
                {
                    regret_ext = regret_array(t).row(d);
                    regret_ext = regret_ext.t();
                    regret_tmp.row(d) = regret_ext.t();
                }
            }

#pragma omp parallel for
            for (unsigned int k = 0; k < K; k++)
            {
                regret.slice(k) = basis_mv(params["basis_mv_idx"](x) - 1).t() * regret_tmp.slice(k) * basis_pr(params["basis_pr_idx"](x) - 1);
                regret.slice(k) *= double(basis_pr(params["basis_pr_idx"](x) - 1).n_cols) / double(P);
                regret.slice(k) *= double(basis_mv(params["basis_mv_idx"](x) - 1).n_cols) / double(D);
            }

            clock.tock("regret");
            clock.tick("learning");
#pragma omp parallel for collapse(2)
            for (unsigned int dr = 0; dr < regret.n_rows; dr++)
            {
                for (unsigned int pr = 0; pr < regret.n_cols; pr++)
                {

// Some aliases for convenience
#define eta_ eta(x).tube(dr, pr)
#define V_ V(x).tube(dr, pr)
#define E_ E(x).tube(dr, pr)
#define r regret.tube(dr, pr)
#define R_ R(x).tube(dr, pr)
#define beta_ beta(x).tube(dr, pr)
#define forget params["forget_regret"](x)
#define beta0 beta0field(x).tube(dr, pr)

                    // Check if a valid method was chosen
                    if (method != "boa" &&
                        method != "bewa" &&
                        method != "ewa" &&
                        method != "ml_poly" &&
                        method != "potential")

                    {
                        Rcpp::stop("Choose 'boa', 'bewa', 'ml_poly' or 'ewa' as method.");
                    }

                    // Learning Rate
                    if (method == "potential")
                    {
                        // TODO
                        eta_ = 1 / (1 / vr(eta_) + square(vr(r)));
                    }
                    else if (method == "ewa")
                    {
                        eta_.fill(params["gamma"](x));
                    }
                    else if (method == "ml_poly")
                    {
                        eta_ = 1 / (1 / vr(eta_) + square(vr(r)));
                    }
                    else if (method == "boa" || method == "bewa")
                    {
                        V_ = vr(V_) * (1 - forget) + square(vr(r));

                        E_ = max(vr(E_) * (1 - forget), abs(vr(r)));

                        eta_ =
                            pmin_arma(
                                min(1 / (2 * vr(E_)),
                                    sqrt(-log(vr(beta0)) / vr(V_))),
                                exp(350));
                    }

                    // Regret
                    if (method == "potential")
                    {
                        R_ = vr(R_ * (1 - forget)) + vr(r);
                    }
                    else if (method == "ewa")
                    {
                        R_ = vr(R_ * (1 - forget)) + vr(r);
                    }
                    else if (method == "ml_poly")
                    {
                        R_ = vr(R_ * (1 - forget)) + vr(r);
                    }
                    else if (method == "boa" || method == "bewa")
                    {

                        vec r_reg = vr(r) - vr(eta_) % square(vr(r));

                        R_ *= (1 - forget); // forget
                        R_ +=
                            0.5 * (r_reg + vr(vr(eta_) % vr(r) > 0.5) % (2 * vr(E_)));
                    }

                    // Weights
                    if (method == "potential")
                    {
                        // TODO
                        beta_ = vr(eta_).t() % vr(R_).t();
                        // beta_ = vr(R_).t();

                        // beta_ /= accu(beta_);
                    }
                    else if (method == "ewa")
                    {
                        beta_ = vr(beta0).t() * K % softmax_r(params["gamma"](x) * vr(R_).t());
                    }
                    else if (method == "ml_poly")
                    {

                        beta_ = vr(beta0).t() * K * params["gamma"](x) % vr(eta_).t() % pmax_arma(vr(R_).t(), exp(-700));
                        beta_ /= accu(beta_);
                    }
                    else if (method == "boa")
                    {
                        // Wintenberger
                        beta_ = vr(beta0).t() * K % softmax_r(log(params["gamma"](x) * vr(eta_).t()) + params["gamma"](x) * vr(eta_).t() % vr(R_).t());
                    }
                    else if (method == "bewa")
                    {
                        // Gaillard
                        beta_ = vr(beta0).t() * K % softmax_r(params["gamma"](x) * vr(eta_).t() % vr(R_).t());
                    }

                    //   // Apply thresholds
                    if (params["soft_threshold"](x) > 0)
                    {
                        int best_k = beta_.index_max();

                        for (double &e : beta_)
                        {
                            threshold_soft(e, params["soft_threshold"](x));
                        }
                        if (accu(beta_) == 0)
                        {
                            beta(x)(dr, pr, best_k) = 1;
                        }
                    }

                    if (params["hard_threshold"](x) > 0)
                    {
                        int best_k = beta_.index_max();
                        for (double &e : beta_)
                        {
                            threshold_hard(e, params["hard_threshold"](x));
                        }
                        if (accu(beta_) == 0)
                        {
                            beta(x)(dr, pr, best_k) = 1;
                        }
                    }

                    // Add fixed_share
                    beta_ =
                        (1 - params["fixed_share"](x)) * vr(beta_) +
                        (params["fixed_share"](x) / K);
                } // pr
            }     // dr
            clock.tock("learning");

#pragma omp parallel for
            // Smoothing
            for (unsigned int k = 0; k < K; k++)
            {
                clock.tick("smoothing");
                weights_tmp(x).slice(k) = hat_mv(params["hat_mv_idx"](x) - 1) *
                                          basis_mv(params["basis_mv_idx"](x) - 1) *
                                          beta(x).slice(k) *
                                          basis_pr(params["basis_pr_idx"](x) - 1).t() *
                                          hat_pr(params["hat_pr_idx"](x) - 1);
                clock.tock("smoothing");
            }

            //  Enshure that constraints hold
            if (method != "potential")
            {
#pragma omp parallel for
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
            }
            if (save_past_performance)
                past_performance(t).slice(x) = tmp_performance;
            R_CheckUserInterrupt();

            // Apply forget
            cum_performance(x) *= (1 - forget_past_performance);
            // Add new loss
            cum_performance(x) += accu(tmp_performance) / double(D * P);

        } // x

        opt_index(t + 1) = cum_performance.index_min();
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
            vec tmp_preds_vec = sum(weights_temp % experts_mat, 1);

            // Sort predictions if quantile_crossing is prohibited
            if (!allow_quantile_crossing)
            {
                tmp_preds_vec = arma::sort(tmp_preds_vec, "ascend", 0);
            }
            predictions.tube(t, d) = tmp_preds_vec;
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
    clock.tock("core");
}

Rcpp::List conline::output()
{
    clock.tick("wrangle");

    // 1-Indexing for R-Output
    opt_index += 1;

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
        Rcpp::Named("allow_quantile_crossing") = allow_quantile_crossing,
        Rcpp::Named("save_past_performance") = save_past_performance,
        Rcpp::Named("save_predictions_grid") = save_predictions_grid);

    Rcpp::List model_objects = Rcpp::List::create(
        Rcpp::Named("weights_tmp") = weights_tmp,
        Rcpp::Named("predictions_grid") = predictions_grid,
        Rcpp::Named("cum_performance") = cum_performance,
        Rcpp::Named("hat_pr") = hat_pr,
        Rcpp::Named("hat_mv") = hat_mv,
        Rcpp::Named("basis_pr") = basis_pr,
        Rcpp::Named("basis_mv") = basis_mv,
        Rcpp::Named("V") = V,
        Rcpp::Named("E") = E,
        Rcpp::Named("eta") = eta,
        Rcpp::Named("R") = R,
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
        Rcpp::Named("opt_index") = opt_index,
        Rcpp::Named("parametergrid") = params,
        Rcpp::Named("params_basis_pr") = params_basis_pr,
        Rcpp::Named("params_basis_mv") = params_basis_mv,
        Rcpp::Named("params_hat_pr") = params_hat_pr,
        Rcpp::Named("params_hat_mv") = params_hat_mv,
        Rcpp::Named("specification") = model_spec);

    out.attr("class") = "online";

    clock.tock("wrangle");
    // Rcpp::List out;
    return out;
}

void conline::init_update(
    Rcpp::List &object,
    arma::mat &new_y,
    arma::field<arma::cube> &new_experts)
{
    clock.tick("init");

    // This creates a references not copies
    Rcpp::List specification = object["specification"];
    Rcpp::List model_parameters = specification["parameters"];
    Rcpp::List model_data = specification["data"];
    Rcpp::List model_objects = specification["objects"];

    // Data

    // Join old and new expert_predictions
    arma::field<cube> old_experts = model_data["experts"];
    experts.set_size(old_experts.n_rows + new_experts.n_rows);
    experts.rows(0, old_experts.n_rows - 1) = old_experts;
    if (new_experts.n_rows > 0)
    {
        experts.rows(old_experts.n_rows, experts.n_rows - 1) = new_experts;
    }

    y = Rcpp::as<arma::mat>(model_data["y"]);
    y.insert_rows(y.n_rows, new_y);
    start = T - new_y.n_rows;

    if (T_E_Y < 0)
        Rcpp::stop("Number of provided expert predictions has to match or exceed observations.");

    tau = Rcpp::as<arma::vec>(model_data["tau"]);

    params = Rcpp::as<std::map<std::string, arma::colvec>>(object["parametergrid"]);
    params_basis_pr = Rcpp::as<std::map<std::string, arma::colvec>>(object["params_basis_pr"]);
    params_basis_mv = Rcpp::as<std::map<std::string, arma::colvec>>(object["params_basis_mv"]);
    params_hat_pr = Rcpp::as<std::map<std::string, arma::colvec>>(object["params_hat_pr"]);
    params_hat_mv = Rcpp::as<std::map<std::string, arma::colvec>>(object["params_hat_mv"]);

    opt_index = Rcpp::as<arma::vec>(object["opt_index"]);
    // Zero indexing in C++
    opt_index -= 1;
    opt_index.resize(T + 1);

    tmp_performance.zeros(D, P);
    cum_performance = Rcpp::as<arma::vec>(model_objects["cum_performance"]);

    weights_tmp =
        Rcpp::as<arma::field<cube>>(model_objects["weights_tmp"]);

    // // Output Objects
    predictions = Rcpp::as<arma::cube>(object["predictions"]);
    predictions.resize(T + T_E_Y, D, P);
    weights.set_size(T + 1);
    weights.rows(0, start) = Rcpp::as<arma::field<cube>>(object["weights"]);

    basis_pr = Rcpp::as<arma::field<arma::sp_mat>>(model_objects["basis_pr"]);
    basis_mv = Rcpp::as<arma::field<arma::sp_mat>>(model_objects["basis_mv"]);
    hat_pr = Rcpp::as<arma::field<arma::sp_mat>>(model_objects["hat_pr"]);
    hat_mv = Rcpp::as<arma::field<arma::sp_mat>>(model_objects["hat_mv"]);

    V = Rcpp::as<arma::field<cube>>(model_objects["V"]);
    E = Rcpp::as<arma::field<cube>>(model_objects["E"]);
    eta = Rcpp::as<arma::field<cube>>(model_objects["eta"]);
    R = Rcpp::as<arma::field<cube>>(model_objects["R"]);
    beta = Rcpp::as<arma::field<cube>>(model_objects["beta"]);
    beta0field = Rcpp::as<arma::field<cube>>(model_objects["beta0field"]);

    // //   // Misc parameters
    lead_time = model_parameters["lead_time"];
    loss_function = Rcpp::as<std::string>(model_parameters["loss_function"]);
    loss_parameter = model_parameters["loss_parameter"];
    loss_gradient = model_parameters["loss_gradient"];
    method = Rcpp::as<std::string>(model_parameters["method"]);

    forget_past_performance = model_parameters["forget_past_performance"];
    allow_quantile_crossing = model_parameters["allow_quantile_crossing"];

    save_past_performance = model_parameters["save_past_performance"];
    save_predictions_grid = model_parameters["save_predictions_grid"];

    if (save_past_performance)
    {
        past_performance.set_size(T);
        past_performance.rows(0, start - 1) =
            Rcpp::as<arma::field<cube>>(object["past_performance"]);
    }
    else
    {
        past_performance.set_size(1);
    }

    if (save_predictions_grid)
    {
        predictions_grid.set_size(T + T_E_Y);
        predictions_grid.rows(0, old_experts.n_rows - 1) = Rcpp::as<arma::field<cube>>(model_objects["predictions_grid"]);
    }
    else
    {
        predictions_grid.set_size(1 + lead_time);
        predictions_grid = Rcpp::as<arma::field<cube>>(model_objects["predictions_grid"]);
    }

    loss_for.zeros(T, D, P);
    loss_for.rows(0, start - 1) =
        Rcpp::as<arma::cube>(object["forecaster_loss"]);

    loss_exp.set_size(T);
    loss_exp.rows(0, start - 1) =
        Rcpp::as<arma::field<cube>>(object["experts_loss"]);

    clock.tock("init");
}
