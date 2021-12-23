//  Student.cpp
#include "conline.h"

// Constructor

// Getters

// Methods

// Functions
// conline class was exposed via "profoc_types.h"
// So we can use it here as input and output if necessary

// [[Rcpp::export]]
conline makeOnline()
{
    conline d;
    return d;
}

// [[Rcpp::export]]
void addtox(conline &x, double &z)
{
    x.x += z;
}