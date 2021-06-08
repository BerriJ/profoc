// [[Rcpp::depends(RcppArmadillo)]]
#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
#define OPTIM_ENABLE_ARMA_WRAPPERS
#include <optim.hpp>
using namespace arma;

using optim::Mat_t;
using optim::Vec_t;

#include <loss.h>

// additional data for objective function
struct objective_data
{
    int K;
    vec truth;
    mat experts;
    double tau;
    std::string loss_function;
    double loss_scaling;
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
    }

    obj_val = mean(loss_vec);

    if (grad_out)
    {
        for (int k = 0; k < K; k++)
        {
            (*grad_out)(k) = loss_grad_wrt_w(experts.col(k),
                                             pred,
                                             truth,
                                             tau,
                                             loss_function,
                                             loss_scaling,
                                             vals_inp(k));
        }
    }

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
vec optimize_weights(vec initvals,
                     const vec &truth,
                     const mat &experts,
                     const bool &convex_constraint = false,
                     const std::string &loss_function = "quantile",
                     const double &tau = 0.5,
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

    bool success;
    optim::algo_settings_t settings;

    if (convex_constraint)
    {
        constraint_data opt_constr_data;

        Mat_t convexity(2, K, fill::ones);
        convexity.row(1) = -convexity.row(1);
        Mat_t constr_positivity(K, K, fill::eye);
        opt_constr_data.constr_mat = join_cols(convexity, -constr_positivity);

        Vec_t constr_rhs(K + 2, fill::zeros);
        constr_rhs(0) = 1.0;
        constr_rhs(1) = -1.0;
        opt_constr_data.constr_rhs = std::move(constr_rhs);

        success = optim::sumt(initvals, objective, &opt_obj_data, constraint_function, &opt_constr_data, settings);
    }
    else
    {
        settings.rel_objfn_change_tol = 1E-07;
        success = optim::nm(initvals, objective, &opt_obj_data, settings);
    }

    if (!success)
    {
        Rcpp::warning("Warning: Convergence was not succesfull.");
    }

    return initvals;
}

// [[Rcpp::export]]
Rcpp::List oracle(mat &y,
                  const cube &experts,
                  Rcpp::NumericVector tau = Rcpp::NumericVector::create(),
                  const std::string loss_function = "quantile",
                  const double &loss_parameter = 1,
                  const bool &convex_constraint = false)
{
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

    // Iinit weights
    vec init_weights(K, fill::ones);
    init_weights /= accu(init_weights);

    // Matrix to store oracle weights
    mat weights(P, K);
    vec oracles_loss(P);

    mat predictions(T, P);

    for (int p = 0; p < P; p++)
    {

        weights.row(p) = optimize_weights(init_weights,
                                          y.col(p),
                                          experts.col(p),
                                          convex_constraint,
                                          loss_function,
                                          tau_vec(p),
                                          loss_parameter)
                             .t();

        if (convex_constraint)
            init_weights = weights.row(p).t();

        mat exp_tmp = experts.col(p);
        predictions.col(p) = exp_tmp * weights.row(p).t();
        oracles_loss(p) = obj_val;
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