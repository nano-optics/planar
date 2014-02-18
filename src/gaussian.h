#ifndef _gaussian_GAUSSIAN_H
#define _gaussian_GAUSSIAN_H

#include <RcppArmadillo.h>

void progress_bar(double x, double N);
arma::mat rotation_y(const double alpha);
arma::mat rotation_z(const double delta);
arma::cx_colvec incident_field(const double psi);
arma::cx_colvec incident_field2(const double psi, const arma::colvec& s1);


#endif
