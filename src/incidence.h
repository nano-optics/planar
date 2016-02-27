#ifndef _incidence_INCIDENCE_H
#define _incidence_INCIDENCE_H

#include <RcppArmadillo.h>

arma::cx_colvec incident_field(const double psi);
arma::cx_colvec incident_field2(const double psi, const arma::colvec& s1);


#endif
