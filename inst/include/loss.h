#ifndef loss_h
#define loss_h

#include <string>

double loss(const double &y,
            const double &x,
            const double &pred,
            const std::string method,
            const double &tau,
            const double &a,
            const bool &gradient);

double loss_grad_wrt_w(const double &expert,
                       const double &pred,
                       const double &truth,
                       const double &tau,
                       const std::string &loss_function,
                       const double &a,
                       const double &w);

#endif
