// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include <iostream>

#include "utils.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma ;
using namespace std;

void progress_bar(double x, double N)
  {
    // how wide you want the progress meter to be
    int totaldotz=40;
    double fraction = x / N;
    // part of the progressmeter that's already "full"
    int dotz = round(fraction * totaldotz);

    // create the "meter"
    int ii=0;
     Rprintf("%3.0f%% [",fraction*100);
    // part  that's full already
    for ( ; ii < dotz;ii++) {
      Rcpp::Rcout << "=";
    }
    // remaining part (spaces)
    for ( ; ii < totaldotz;ii++) {
       Rprintf(" ");
    }
    // and back to line begin -
    // do not forget the fflush to avoid output buffering problems!
    Rprintf("]\r");
    // fflush(stdout);
  }

arma::mat rotation_y(const double alpha)
  {
    arma::mat Rot(3,3);
    const double ca = cos(alpha), sa = sin(alpha);

    Rot(0,0) = ca;
    Rot(0,1) = 0;
    Rot(0,2) = sa;

    Rot(1,0) = 0;
    Rot(1,1) = 1;
    Rot(1,2) = 0;

    Rot(2,0) = -sa;
    Rot(2,1) = 0;
    Rot(2,2) = ca;
    return (Rot);
  }

arma::mat rotation_z(const double delta)
  {
    arma::mat Rot(3,3);
    const double ca = cos(delta), sa = sin(delta);

    Rot(0,0) = ca;
    Rot(0,1) = sa;
    Rot(0,2) = 0;

    Rot(1,0) = -sa;
    Rot(1,1) = ca;
    Rot(1,2) = 0;

    Rot(2,0) = 0;
    Rot(2,1) = 0;
    Rot(2,2) = 1;
    return (Rot);
  }
