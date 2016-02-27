#ifndef _utils_UTILS_H
#define _utils_UTILS_H

#include <RcppArmadillo.h>

void progress_bar(double x, double N);
arma::mat rotation_y(const double alpha);
arma::mat rotation_z(const double delta);


#endif
