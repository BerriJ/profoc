// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <string>
#include <progress.hpp>
#include <loss.h>

using namespace arma;

// [[Rcpp::export]]
double loss(const double &y,
            const double &x,
            const double &pred = 0,
            const std::string method = "quantile",
            const double &tau = 0.5,
            const double &a = 1,
            const bool &gradient = true)
{
    double loss;

    if (method == "quantile")
    {
        if (!gradient)
        {
            loss = ((y < x) - tau) * (sign(x) * std::pow(std::abs(x), a) - sign(y) * std::pow(std::abs(y), a));
        }
        else
        {
            loss = ((pred >= y) - tau) * (a * std::pow(std::abs(pred), (a - 1)));
            loss = loss * x;
        }
    }
    else if (method == "expectile")
    {
        if (!gradient)
        {
            loss = 2 * std::abs((x >= y) - tau) * (std::pow(std::abs(y), (a + 1)) - std::pow(std::abs(x), (a + 1)) - (a + 1) * sign(x) * std::pow(std::abs(x), a) * (y - x));
        }
        else
        {
            loss = 2 * std::abs((pred >= y) - tau) * (-a * (a + 1) * (y - pred) * std::pow(std::abs(pred), (a - 1)));
            loss = loss * x;
        }
    }
    else if (method == "percentage")
    {
        if (!gradient)
        {
            loss = std::abs(1 - std::pow(x / y, a));
        }
        else
        {
            loss = a * (std::pow(pred / y, a) - 1) * std::pow(pred / y, a) / (pred * std::abs(1 - std::pow(pred / y, a)));
            loss = loss * x;
        }
    }
    else
    {
        Rcpp::stop("Choose quantile loss 'quantile' expectiles 'expectile' or as 'percentage' loss.");
    }
    return loss;
}

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

// Loss gradient w.r.t. weights
// [[Rcpp::export]]
double loss_grad_wrt_w(const double &expert,
                       const double &pred,
                       const double &truth,
                       const double &tau,
                       const std::string &loss_function,
                       const double &a,
                       const double &w)
{

    double loss_grad;

    if (loss_function == "quantile")
    {
        loss_grad = a * expert * pow(fabs(pred), a - 1) * ((pred >= truth) - tau);
    }
    else if (loss_function == "expectile")
    {
        loss_grad = sgn(pred - truth) * expert * (pred - truth) * (-2 * tau + 2 * (pred >= truth));
    }
    else if (loss_function == "percentage")
    {
        double nom = a * w * pow(pred / truth, a - 1) * (1 - pow(pred / truth, a));
        double denom = truth * abs(1 - pow(pred / truth, a));

        loss_grad = -nom / denom;
    }
    else
    {
        Rcpp::stop("Choose quantile loss 'quantile' expectiles 'expectile' or as 'percentage' loss.");
    }

    return loss_grad;
}