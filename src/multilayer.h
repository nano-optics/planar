#ifndef _multilayer_MULTILAYER_H
#define _multilayer_MULTILAYER_H

#include <RcppArmadillo.h>

Rcpp::List cpp_multilayer_field(const double k0,				
			    const double kx,
			    const arma::cx_vec& epsilon,
			    const arma::colvec& thickness,
			    const arma::colvec& z, const double psi);

#endif
