#include <misc.h>
#include <loss.h>

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <clock.h>
#include <thread>

//  Student.cpp
#include "conline.h"

// conline class was exposed via "profoc_types.h"
// So we can use it here as input and output if necessary

// Constructor

// Getters

// Methods

// Functions
void conline::init_objects()
{
    // Expand tau if necessary
    if (tau.n_elem == 1)
    {
        tau.resize(P);
        tau.fill(tau(0));
    }
}

// // [[Rcpp::export]]
// arma::field<arma::mat> test()
// {
//     arma::mat A(3, 3);
//     arma::field<arma::mat> F(1, 1, 1);
//     F(0, 0, 0) = A;
//     return (F);
// }