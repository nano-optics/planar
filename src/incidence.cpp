// [[Rcpp::depends(RcppArmadillo)]]



#include <RcppArmadillo.h>
#include <iostream>

#include "multilayer.h"
#include "cubature.h"
#include "utils.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma ;
using namespace std;


arma::cx_colvec incident_field(const double psi)
  {
   arma::cx_colvec E(3);

    E(0) = cos(psi);
    E(1) = sin(psi);
    E(2) = 0;

    return (E);
  }


arma::cx_colvec incident_field2(const double psi, const arma::colvec& s1)
  {
   arma::cx_colvec E(3);
   double s1z = s1(2), s1x=s1(0), s1y=s1(1);
   double s1z2 = s1z*s1z, s1y2 = s1y*s1y;

   E(0) = cos(psi)*(s1z + s1y2*(1-s1z)/(1-s1z2)) + sin(psi)*((s1z - 1)*s1x*s1y/(1-s1z2));
   E(1) = cos(psi)*((s1z - 1)*s1x*s1y/(1-s1z2)) + sin(psi)*(1-s1y2*(1-s1z)/(1-s1z2));
   E(2) = -s1x*cos(psi) -s1y*sin(psi);

   return (E);
  }
