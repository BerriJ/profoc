#include <RcppArmadillo.h>
#include <progress.hpp>
#include <thread>

#include "misc.h"
#include "loss.h"
#include "rcpptimer.h"
#include "conline.h"
#include "profoc_types.h"

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
  timer.autoreturn = get_timings;

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

  predictions_got_sorted.zeros(T + T_E_Y, D);
}

void conline::set_grid_objects()
{
  timer.tic("init");

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

  for (unsigned int x = 0; x < X; x++)
  {
    unsigned int Pr = basis_pr(params["basis_pr_idx"](x) - 1).n_cols;
    unsigned int Dr = basis_mv(params["basis_mv_idx"](x) - 1).n_cols;

    // Learning parameters
    V(x).zeros(Dr, Pr, K);
    E(x).zeros(Dr, Pr, K);

    arma::cube eta_(Dr, Pr, K, arma::fill::zeros);
    eta(x) = eta_;
    if (method == "ml_poly")
    {
      eta_.fill(std::exp(350));
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
      beta(x).slice(k) = arma::pinv(
                             arma::mat(basis_mv(params["basis_mv_idx"](x) - 1))) *
                         w0.slice(k) *
                         arma::pinv(arma::mat(basis_pr(params["basis_pr_idx"](x) - 1))).t();
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
      past_performance(t).fill(arma::datum::nan);
    }
    predictions_grid(t).set_size(D, P, X);

    // Store predictions w.r.t. grid for time t
    arma::cube tmp_preds_cube(D, P, X);

    for (unsigned int d = 0; d < D; d++)
    {
      // Save final weights weights_tmp
      weights(t).row(d) = weights_tmp(opt_index(t)).row(d);

      // Store expert predictions temporarily
      arma::mat experts_mat = experts(t).row(d);

      for (unsigned int x = 0; x < X; x++)
      {

        arma::mat weights_temp = weights_tmp(x).row(d);

        // Forecasters prediction
        arma::vec tmp_preds_vec = arma::sum(weights_temp % experts_mat, 1);

        tmp_preds_cube(arma::span(d), arma::span::all, arma::span(x)) = tmp_preds_vec;
      }
    }
    predictions_grid(t) = tmp_preds_cube;
    // Final Prediction
    predictions.row(t) = tmp_preds_cube.slice(opt_index(t));
  }

  start = lead_time;

  timer.toc("init");
}

void conline::learn()
{
  Progress prog(T, trace);

  timer.tic("core");

  for (unsigned int tp = 0; tp < predictions_grid.n_rows; tp++)
  {
    predictions_grid(tp).set_size(D, P, X);
  }

  for (unsigned int t = start; t < T; t++)
  {

    weights(t).set_size(D, P, K);

    if (save_past_performance)
      past_performance(t).set_size(D, P, X);

    timer.tic("loss");

    // Store predictions w.r.t. grid for time t
    arma::cube tmp_preds_cube(D, P, X);

    for (unsigned int d = 0; d < D; d++)
    {
      // Save final weights weights_tmp
      weights(t).row(d) = weights_tmp(opt_index(t)).row(d);

      // Store expert predictions temporarily
      arma::mat experts_mat = experts(t).row(d);

      // Predictions using different parameter values
      for (unsigned int x = 0; x < X; x++)
      {

        arma::mat weights_temp = weights_tmp(x).row(d);
        arma::vec tmp_preds_vec = arma::sum(weights_temp % experts_mat, 1);

        // Sort predictions if quantile_crossing is prohibited
        if (!allow_quantile_crossing)
        {
          if ((x == opt_index(t)) && (!tmp_preds_vec.is_sorted()))
          {
            predictions_got_sorted(t, d) = 1;
          }
          tmp_preds_vec = arma::sort(tmp_preds_vec, "ascend", 0);
        }
        tmp_preds_cube(arma::span(d), arma::span::all, arma::span(x)) = tmp_preds_vec;
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

    timer.toc("loss");
    for (unsigned int x = 0; x < X; x++)
    {

      timer.tic("regret");
      arma::mat lexp_int(P, K); // Experts loss
      arma::mat lexp_ext(P, K); // Experts loss
      arma::mat lexp(P, K);     // Experts loss
      arma::vec lfor(P);        // Forecasters loss
      arma::cube regret_tmp(D, P, K);
      arma::cube regret(basis_mv(params["basis_mv_idx"](x) - 1).n_cols,
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
              lexp_ext.row(p) = arma::vectorise(loss_array(t).tube(d, p)).t();
              lexp.row(p) = (1 - params["loss_share"](x)) * lexp_int.row(p) + params["loss_share"](x) * lexp_ext.row(p);
            }
          }
          else
          {
            lexp_ext.row(p) = arma::vectorise(loss_array(t).tube(d, p)).t();
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

        arma::mat regret_int(P, K);
        arma::mat regret_ext(P, K);

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

      timer.toc("regret");

      timer.tic("learning");
#pragma omp parallel for collapse(2)
      for (unsigned int dr = 0; dr < regret.n_rows; dr++)
      {
        for (unsigned int pr = 0; pr < regret.n_cols; pr++)
        {

          arma::vec r = regret.tube(dr, pr);

          if (method == "ewa")
          {
            // Update the cumulative regret used by eta
            R(x).tube(dr, pr) = arma::vectorise(R(x).tube(dr, pr) * (1 - params["forget_regret"](x))) + r;
            eta(x).tube(dr, pr).fill(params["gamma"](x));
            beta(x).tube(dr, pr) = arma::vectorise(beta0field(x).tube(dr, pr)).t() * K % softmax_r(params["gamma"](x) * arma::vectorise(R(x).tube(dr, pr)).t());
          }
          else if (method == "ml_poly")
          {
            // Update the cumulative regret used by ML_Poly
            R(x)
                .tube(dr, pr) = arma::vectorise(R(x).tube(dr, pr) * (1 - params["forget_regret"](x))) + r;

            // Update the learning rate
            eta(x).tube(dr, pr) = 1 / (1 / arma::vectorise(eta(x).tube(dr, pr)).t() + arma::square(r.t()));

            beta(x).tube(dr, pr) = arma::vectorise(beta0field(x).tube(dr, pr)).t() * K * params["gamma"](x) % arma::vectorise(eta(x).tube(dr, pr)).t() % pmax_arma(arma::vectorise(R(x).tube(dr, pr)).t(), std::exp(-700));
            beta(x).tube(dr, pr) /= arma::accu(beta(x).tube(dr, pr));
          }
          else if (method == "boa" || method == "bewa")
          {
            V(x).tube(dr, pr) = arma::vectorise(V(x).tube(dr, pr)).t() * (1 - params["forget_regret"](x)) + arma::square(r.t());

            E(x).tube(dr, pr) = pmax_arma(arma::max(arma::vectorise(E(x).tube(dr, pr)).t() * (1 - params["forget_regret"](x)), arma::abs(r.t())), std::exp(-350));

            eta(x)
                .tube(dr, pr) =
                pmin_arma(
                    arma::min(1 / (2 * arma::vectorise(E(x).tube(dr, pr))),
                              arma::sqrt(-arma::log(arma::vectorise(beta0field(x).tube(dr, pr))) / pmax_arma(arma::vectorise(V(x).tube(dr, pr)), std::exp(-350)))),
                    std::exp(350));

            arma::vec r_reg = r - arma::vectorise(eta(x).tube(dr, pr)) % arma::square(r);

            R(x).tube(dr, pr) *= (1 - params["forget_regret"](x)); // forget
            R(x).tube(dr, pr) +=
                0.5 * (r_reg + arma::conv_to<arma::colvec>::from(arma::vectorise(eta(x).tube(dr, pr)) % r > 0.5) % (2 * arma::vectorise(E(x).tube(dr, pr))));

            if (method == "boa")
            {
              // Wintenberger
              beta(x).tube(dr, pr) = arma::vectorise(beta0field(x).tube(dr, pr)).t() * K % softmax_r(arma::log(params["gamma"](x) * arma::vectorise(eta(x).tube(dr, pr)).t()) + params["gamma"](x) * arma::vectorise(eta(x).tube(dr, pr)).t() % arma::vectorise(R(x).tube(dr, pr)).t());
            }
            else
            {
              // Gaillard
              beta(x).tube(dr, pr) = arma::vectorise(beta0field(x).tube(dr, pr)).t() * K % softmax_r(params["gamma"](x) * arma::vectorise(eta(x).tube(dr, pr)).t() % arma::vectorise(R(x).tube(dr, pr)).t());
            }
          }
          else
          {
            Rcpp::stop("Choose 'boa', 'bewa', 'ml_poly' or 'ewa' as method.");
          }

          //   // Apply thresholds
          if (params["soft_threshold"](x) > 0)
          {
            int best_k = beta(x).tube(dr, pr).index_max();

            for (double &e : beta(x).tube(dr, pr))
            {
              threshold_soft(e, params["soft_threshold"](x));
            }
            if (arma::accu(beta(x).tube(dr, pr)) == 0)
            {
              beta(x)(dr, pr, best_k) = 1;
            }
          }

          if (params["hard_threshold"](x) > 0)
          {
            int best_k = beta(x).tube(dr, pr).index_max();
            for (double &e : beta(x).tube(dr, pr))
            {
              threshold_hard(e, params["hard_threshold"](x));
            }
            if (arma::accu(beta(x).tube(dr, pr)) == 0)
            {
              beta(x)(dr, pr, best_k) = 1;
            }
          }

          // Add fixed_share
          beta(x).tube(dr, pr) =
              (1 - params["fixed_share"](x)) * arma::vectorise(beta(x).tube(dr, pr)) +
              (params["fixed_share"](x) / K);
        } // pr
      } // dr

      timer.toc("learning");

#pragma omp parallel for
      // Smoothing
      for (unsigned int k = 0; k < K; k++)
      {

        timer.tic("smoothing");
        weights_tmp(x).slice(k) = hat_mv(params["hat_mv_idx"](x) - 1) *
                                  basis_mv(params["basis_mv_idx"](x) - 1) *
                                  beta(x).slice(k) *
                                  basis_pr(params["basis_pr_idx"](x) - 1).t() *
                                  hat_pr(params["hat_pr_idx"](x) - 1);

        timer.toc("smoothing");
      }

#pragma omp parallel for
      //  Enshure that constraints hold
      for (unsigned int p = 0; p < P; p++)
      {
        for (unsigned int d = 0; d < D; d++)
        {
          // Positivity
          weights_tmp(x)(arma::span(d), arma::span(p), arma::span::all) =
              pmax_arma(weights_tmp(x)(arma::span(d), arma::span(p), arma::span::all), std::exp(-700));

          // // Affinity
          weights_tmp(x)(arma::span(d), arma::span(p), arma::span::all) /=
              arma::accu(weights_tmp(x)(arma::span(d), arma::span(p), arma::span::all));
        }
      }
      if (save_past_performance)
        past_performance(t).slice(x) = tmp_performance;
      R_CheckUserInterrupt();

      // Apply forget
      cum_performance(x) *= (1 - forget_past_performance);
      // Add new loss
      cum_performance(x) += arma::accu(tmp_performance) / double(D * P);

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
      arma::mat experts_mat = experts(t).row(d);
      arma::mat weights_temp = weights(T).row(d);
      arma::vec tmp_preds_vec = arma::sum(weights_temp % experts_mat, 1);

      // Sort predictions if quantile_crossing is prohibited
      if (!allow_quantile_crossing)
      {
        if (!tmp_preds_vec.is_sorted())
        {
          predictions_got_sorted(t, d) = 1;
        }
        tmp_preds_vec = arma::sort(tmp_preds_vec, "ascend", 0);
      }
      predictions.tube(t, d) = tmp_preds_vec;
    }
  }

  // Save losses suffered by forecaster and experts

  timer.tic("loss_for_exp");
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

  timer.toc("loss_for_exp");

  timer.toc("core");
}

void conline::init_update(
    Rcpp::List &object,
    arma::mat &new_y,
    arma::field<arma::cube> &new_experts)
{

  timer.tic("init update");

  // This creates a references not copies
  Rcpp::List specification = object["specification"];
  Rcpp::List model_parameters = specification["parameters"];
  Rcpp::List model_data = specification["data"];
  Rcpp::List model_objects = specification["objects"];

  // Data

  // Join old and new expert_predictions
  arma::field<arma::cube> old_experts = model_data["experts"];
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
  {
    timer.toc("init update");
    Rcpp::stop("Number of provided expert predictions has to match or exceed observations.");
  }

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
      Rcpp::as<arma::field<arma::cube>>(model_objects["weights_tmp"]);

  // // Output Objects
  predictions = Rcpp::as<arma::cube>(object["predictions"]);
  predictions.resize(T + T_E_Y, D, P);
  predictions_got_sorted = Rcpp::as<arma::mat>(object["predictions_got_sorted"]);
  predictions_got_sorted.resize(T + T_E_Y, D);
  weights.set_size(T + 1);
  weights.rows(0, start) = Rcpp::as<arma::field<arma::cube>>(object["weights"]);

  basis_pr = Rcpp::as<arma::field<arma::sp_mat>>(model_objects["basis_pr"]);
  basis_mv = Rcpp::as<arma::field<arma::sp_mat>>(model_objects["basis_mv"]);
  hat_pr = Rcpp::as<arma::field<arma::sp_mat>>(model_objects["hat_pr"]);
  hat_mv = Rcpp::as<arma::field<arma::sp_mat>>(model_objects["hat_mv"]);

  V = Rcpp::as<arma::field<arma::cube>>(model_objects["V"]);
  E = Rcpp::as<arma::field<arma::cube>>(model_objects["E"]);
  eta = Rcpp::as<arma::field<arma::cube>>(model_objects["eta"]);
  R = Rcpp::as<arma::field<arma::cube>>(model_objects["R"]);
  beta = Rcpp::as<arma::field<arma::cube>>(model_objects["beta"]);
  beta0field = Rcpp::as<arma::field<arma::cube>>(model_objects["beta0field"]);

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
        Rcpp::as<arma::field<arma::cube>>(object["past_performance"]);
  }
  else
  {
    past_performance.set_size(1);
  }

  if (save_predictions_grid)
  {
    predictions_grid.set_size(T + T_E_Y);
    predictions_grid.rows(0, old_experts.n_rows - 1) = Rcpp::as<arma::field<arma::cube>>(model_objects["predictions_grid"]);
  }
  else
  {
    predictions_grid.set_size(1 + lead_time);
    predictions_grid = Rcpp::as<arma::field<arma::cube>>(model_objects["predictions_grid"]);
  }

  loss_for.zeros(T, D, P);
  loss_for.rows(0, start - 1) =
      Rcpp::as<arma::cube>(object["forecaster_loss"]);

  loss_exp.set_size(T);
  loss_exp.rows(0, start - 1) =
      Rcpp::as<arma::field<arma::cube>>(object["experts_loss"]);

  timer.toc("init update");
}
