#ifndef loss_h
#define loss_h

#include <common_header.h>

double loss(const double &y,
            const double &x,
            const double &pred = 0,
            const std::string method = "quantile",
            const double &tau = 0.5,
            const double &a = 1,
            const bool &gradient = true);

#endif