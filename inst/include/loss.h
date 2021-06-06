#ifndef loss_h
#define loss_h

double loss(const double &y,
            const double &x,
            const double &pred,
            const std::string method,
            const double &tau,
            const double &a,
            const bool &gradient);

#endif