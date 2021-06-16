#include <oracle.h>
#include <misc.h>
#include <loss.h>

// [[Rcpp::depends(RcppArmadillo)]]
#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#define OPTIM_ENABLE_ARMA_WRAPPERS
#include <optim.hpp>

#include <RcppArmadillo.h>

using namespace arma;

using optim::Mat_t;
using optim::Vec_t;

// additional data for objective function
struct objective_data
{
    int K;
    vec truth;
    mat experts;
    double tau;
    std::string loss_function;
    double loss_scaling;
    double forget;
    double penalty_parameter = 0;
    bool intercept;
    bool debias;
};

// Global decleration of objective val
double obj_val;

double objective(const vec &vals_inp, vec *grad_out, void *opt_data)
{
    objective_data *objfn_data = reinterpret_cast<objective_data *>(opt_data);

    int K = objfn_data->K;
    vec truth = objfn_data->truth;
    mat experts = objfn_data->experts;
    std::string loss_function = objfn_data->loss_function;
    double tau = objfn_data->tau;
    double loss_scaling = objfn_data->loss_scaling;
    double forget = objfn_data->forget;
    double penalty_parameter = objfn_data->penalty_parameter;
    bool intercept = objfn_data->intercept;
    bool debias = objfn_data->debias;

    vec pred = experts * vals_inp;
    vec loss_vec(experts.n_rows);

    for (unsigned int i = 0; i < experts.n_rows; i++)
    {
        loss_vec(i) = loss(truth(i),
                           pred(i),
                           0, // Only used when gradient = TRUE
                           loss_function,
                           tau,
                           loss_scaling, // Scaling parameter
                           false);
        loss_vec(i) *= std::pow(1 - forget, experts.n_rows - (i + 1));
    }

    obj_val = sum(loss_vec);

    double constraint_penalty;

    constraint_penalty = penalty_parameter * fabs(sum(vals_inp.subvec(debias * intercept, vals_inp.n_elem - 1)) - 1);

    vec loss_grad(experts.n_rows);

    if (grad_out)
    {
        for (int k = 0; k < K; k++)
        {
            for (unsigned int i = 0; i < experts.n_rows; i++)
            {
                loss_grad(i) = loss_grad_wrt_w(experts(i, k),
                                               pred(i),
                                               truth(i),
                                               tau,
                                               loss_function,
                                               loss_scaling,
                                               vals_inp(k));
                loss_grad(i) *= std::pow(1 - forget, experts.n_rows - (i + 1));
            }

            (*grad_out)(k) = sum(loss_grad);
        }
    }

    // arma::cout << obj_val << arma::endl;
    // arma::cout << vals_inp << arma::endl;
    // arma::cout << constraint_penalty << arma::endl;

    obj_val += constraint_penalty;

    return obj_val;
}

// additional data for constraint function
struct constraint_data
{
    Mat_t constr_mat;
    Vec_t constr_rhs;
};

Vec_t constraint_function(const Vec_t &vals_inp, Mat_t *jacob_out, void *opt_data)
{

    constraint_data *cnstrn_data = reinterpret_cast<constraint_data *>(opt_data);

    Mat_t constr_mat = cnstrn_data->constr_mat;
    Vec_t constr_rhs = cnstrn_data->constr_rhs;

    Vec_t constr_vals = constr_mat * vals_inp - constr_rhs;

    if (jacob_out)
    {
        OPTIM_MATOPS_SET_SIZE_POINTER(jacob_out, constr_mat.n_rows, constr_mat.n_cols);

        (*jacob_out) = constr_mat;
    }

    return constr_vals;
}

// [[Rcpp::export]]
vec optimize_weights(const vec &truth,
                     const mat &experts,
                     const bool &affine = false,
                     const bool &positive = false,
                     const bool &intercept = false,
                     const bool &debias = true,
                     const std::string &loss_function = "quantile",
                     const double &tau = 0.5,
                     const double &forget = 0,
                     const double &loss_scaling = 1)
{

    const int K = experts.n_cols;

    // Prepare additional data for objective
    objective_data opt_obj_data;
    opt_obj_data.K = K;
    opt_obj_data.truth = std::move(truth);
    opt_obj_data.experts = std::move(experts);
    opt_obj_data.loss_function = std::move(loss_function);
    opt_obj_data.tau = std::move(tau);
    opt_obj_data.loss_scaling = std::move(loss_scaling);
    opt_obj_data.forget = std::move(forget);
    opt_obj_data.intercept = intercept;
    opt_obj_data.debias = debias;

    // Iinit weights
    vec initvals(K, fill::zeros);
    initvals.subvec(debias * intercept, initvals.n_elem - 1).fill(1);
    initvals.subvec(debias * intercept, initvals.n_elem - 1) /= (K - debias * intercept);

    bool success;
    optim::algo_settings_t settings;
    constraint_data opt_constr_data;

    if (affine && positive)
    {
        initvals += 1E-08;
        opt_obj_data.penalty_parameter = 0.1;

        settings.vals_bound = true;
        settings.lower_bounds = OPTIM_MATOPS_ZERO_VEC(K);
        settings.lower_bounds.fill(0);
        if (debias && intercept)
        {
            settings.lower_bounds(0) = -exp(700);
        }
        settings.upper_bounds = OPTIM_MATOPS_ZERO_VEC(K);
        settings.upper_bounds.fill(exp(700));

        while (fabs(sum(initvals.subvec(debias * intercept, initvals.n_elem - 1)) - 1) >= 1E-08)
        {
            if (intercept && debias)
            {
                initvals.subvec(debias * intercept, initvals.n_elem - 1).fill(1);
                initvals.subvec(debias * intercept, initvals.n_elem - 1) /= (K - debias * intercept);
            }
            opt_obj_data.penalty_parameter *= 10;
            success = optim::nm(initvals, objective, &opt_obj_data, settings);
            if (opt_obj_data.penalty_parameter > 1E+06)
            {
                break;
            }
        }
    }
    else if (affine)
    {
        initvals += 1E-08;
        opt_obj_data.penalty_parameter = 0.1;

        while (fabs(sum(initvals.subvec(debias * intercept, initvals.n_elem - 1)) - 1) >= 1E-08)
        {
            if (intercept && debias)
            {
                initvals.subvec(debias * intercept, initvals.n_elem - 1).fill(1);
                initvals.subvec(debias * intercept, initvals.n_elem - 1) /= (K - debias * intercept);
            }
            opt_obj_data.penalty_parameter *= 10;
            success = optim::nm(initvals, objective, &opt_obj_data, settings);
            if (opt_obj_data.penalty_parameter > 1E+06)
            {
                break;
            }
        }
    }
    else if (positive)
    {
        settings.vals_bound = true;
        settings.lower_bounds = OPTIM_MATOPS_ZERO_VEC(K);
        settings.lower_bounds.fill(0);
        if (debias && intercept)
        {
            settings.lower_bounds(0) = -exp(700);
        }
        settings.upper_bounds = OPTIM_MATOPS_ZERO_VEC(K);
        settings.upper_bounds.fill(exp(700));

        success = optim::nm(initvals, objective, &opt_obj_data, settings);
    }
    else
    {
        success = optim::nm(initvals, objective, &opt_obj_data, settings);
    }

    if (!success)
    {
        Rcpp::warning("Warning: Convergence was not succesfull.");
    }

    return initvals;
}

//' @template function_oracle
//'
//' @template param_y
//' @template param_experts
//' @template param_tau
//' @template param_affine
//' @template param_positive
//' @template param_intercept
//' @template param_debias
//' @template param_loss_function
//' @template param_loss_parameter
//' @template param_forget
//' @usage oracle(y, experts, tau, affine = FALSE,
//' positive = FALSE, intercept = FALSE, debias = TRUE,
//' loss_function = "quantile", loss_parameter = 1, forget = 0)
//' @export
// [[Rcpp::export]]
Rcpp::List oracle(arma::mat &y,
                  cube &experts,
                  Rcpp::NumericVector tau = Rcpp::NumericVector::create(),
                  const bool &affine = false,
                  const bool &positive = false,
                  const bool &intercept = false,
                  const bool &debias = true,
                  const std::string loss_function = "quantile",
                  const double &loss_parameter = 1,
                  const double &forget = 0)
{

    if (intercept)
    {
        mat intercept_mat(experts.n_rows, experts.n_cols, fill::ones);
        experts = join_slices(intercept_mat, experts);
    }

    // Object Dimensions
    const int T = y.n_rows;
    const int P = experts.n_cols;
    const int K = experts.n_slices;

    // Preprocessing of inputs
    if (y.n_cols == 1)
    {
        y = arma::repmat(y, 1, P);
    }

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

    // Matrix to store oracle weights
    mat weights(P, K);
    vec oracles_loss(P);

    mat predictions(T, P);

    for (int p = 0; p < P; p++)
    {

        weights.row(p) = optimize_weights(y.col(p),
                                          experts.col(p),
                                          affine,
                                          positive,
                                          intercept,
                                          debias,
                                          loss_function,
                                          tau_vec(p),
                                          forget,
                                          loss_parameter)
                             .t();

        mat exp_tmp = experts.col(p);
        predictions.col(p) = exp_tmp * weights.row(p).t();
        oracles_loss(p) = obj_val / experts.n_rows;
        R_CheckUserInterrupt();
    }

    cube loss_exp(T, P, K, fill::zeros);

    for (int t = 0; t < T; t++)
    {
        for (int p = 0; p < P; p++)
        {
            for (int k = 0; k < K; k++)
            {
                loss_exp(t, p, k) =
                    loss(y(t, p),
                         experts(t, p, k),
                         0,              // where to evaluate the gradient
                         loss_function,  // method
                         tau_vec(p),     // tau_vec
                         loss_parameter, // alpha
                         false);         // gradient
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("predictions") = predictions,
        Rcpp::Named("weights") = weights,
        Rcpp::Named("oracles_loss") = oracles_loss,
        Rcpp::Named("experts_loss") = loss_exp);
}