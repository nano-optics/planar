// [[Rcpp::depends(RcppArmadillo)]]

//
// multilayer reflection and transmission coefficients
// also computes the internal fields
//

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma ;
using namespace std;


// [[Rcpp::export]]
Rcpp::List cpp_layer_fresnel(const arma::colvec& k0,
			     const arma::cx_mat& kx,
			     const arma::cx_mat& epsilon,
			     const double& thickness) {


  const int Nlambda = k0.n_rows;
  const int Ntheta = kx.n_cols;

  const arma::colvec& k02 = k0 % k0;
  const arma::cx_mat& kx2 = kx % kx;

  // kiz
  arma::cx_mat kiz0 = arma::ones<arma::cx_mat>(Nlambda, Ntheta), kiz1=kiz0, kiz2=kiz0;

  kiz0  = sqrt(repmat(epsilon.col(0) % k02, 1, Ntheta)  - kx2);
  kiz1  = sqrt(repmat(epsilon.col(1) % k02, 1, Ntheta)  - kx2);
  kiz2  = sqrt(repmat(epsilon.col(2) % k02, 1, Ntheta)  - kx2);

  // calculate single-interface fresnel coefficients and phase factors

  arma::cx_mat ap = arma::cx_mat(Nlambda, Ntheta), as=ap, bp=ap, bs=ap;
  arma::cx_mat rp01 = arma::cx_mat(Nlambda, Ntheta), tp01=rp01, rp12=rp01, tp12=rp01,\
    rs01=rp01, ts01=rp01, rs12=rp01, ts12=rp01, phase1=rp01, phase2=rp01;

  // p-polarisation

  ap = kiz0 / repmat(epsilon.col(0), 1, Ntheta);
  bp = kiz1 / repmat(epsilon.col(1), 1, Ntheta);
  rp01 = (ap - bp) / (ap + bp);
  tp01 = 2 * ap / (ap + bp);

  ap = kiz1 / repmat(epsilon.col(1), 1, Ntheta);
  bp = kiz2 / repmat(epsilon.col(2), 1, Ntheta);
  rp12 = (ap - bp) / (ap + bp);
  tp12 = 2 * ap / (ap + bp);

   // s-polarisation

  rs01 = (kiz0 - kiz1) / (kiz0 + kiz1);
  ts01 = 2 * kiz0 / (kiz0 + kiz1);
  rs12 = (kiz1 - kiz2) / (kiz1 + kiz2);
  ts12 = 2 * kiz1 / (kiz1 + kiz2);

  phase2 =  exp(arma::cx_double(0,2)*thickness*kiz1);
  phase1 =  exp(arma::cx_double(0,1)*thickness*kiz1);

  arma::cx_mat rp = rp01, rs = rp01, tp = rp01, ts = rp01;

  rp  = ( rp01 + rp12%phase2 ) / ( 1 + rp01%rp12%phase2 );
  rs  = ( rs01 + rs12%phase2 ) / ( 1 + rs01%rs12%phase2 );
  tp  = ( tp01 % tp12%phase1 ) / ( 1 + rp01%rp12%phase2 );
  ts  = ( ts01 % ts12%phase1 ) / ( 1 + rs01%rs12%phase2 );

  return List::create(
		      _["rp"]  = rp,
		      _["rs"]  = rs,
		      _["tp"]  = tp,
		      _["ts"]  = ts
 ) ;

}


// [[Rcpp::export]]
Rcpp::List cpp_recursive_fresnel(const arma::colvec& k0,
			     const arma::cx_mat& kx,
			     const arma::cx_mat& epsilon,
			     const arma::colvec& thickness,
			     const int& polarisation) {


  const int Nlambda = k0.n_rows;
  const int Ntheta = kx.n_cols;
  const int Nlayer = thickness.n_rows;

  const arma::colvec& k02 = k0 % k0;
  const arma::cx_mat& kx2 = kx % kx;

  // loop to calculate kiz
  arma::cx_cube kiz = arma::ones<arma::cx_cube>(Nlambda, Ntheta, Nlayer);

  // loop from first to last medium, define normal kiz
  int ii=0;
  for(ii=0; ii < Nlayer; ii++){
    kiz.slice(ii)  = sqrt(repmat(epsilon.col(ii) % k02, 1, Ntheta)  - kx2);
  }

  // calculate all single-interface fresnel coefficients and phase factors

  const int Nsurface = Nlayer - 1 ;

  arma::cx_cube rsingle = 0*arma::ones<arma::cx_cube>(Nlambda, Ntheta, Nsurface), tsingle=rsingle;
  arma::cx_cube phase1 = arma::ones<arma::cx_cube>(Nlambda, Ntheta, Nlayer), phase2 = phase1;
  arma::cx_mat a = arma::cx_mat(Nlambda, Ntheta), b=a, r=a, t=a;

  // loop over interfaces
  for (ii = 0; ii < Nsurface; ii++){

    if (polarisation == 0){ // p-polarisation

      a = kiz.slice(ii) / repmat(epsilon.col(ii), 1, Ntheta);
      b = kiz.slice(ii + 1) / repmat(epsilon.col(ii + 1), 1, Ntheta);

    } else { // s-polarisation

      a = kiz.slice(ii);
      b = kiz.slice(ii + 1);

    }

    rsingle.slice(ii) = (a - b) / (a + b);
    tsingle.slice(ii) = 2 * a / (a + b) ;

    phase2.slice(ii) =  exp(arma::cx_double(0,2)*thickness(ii)*kiz.slice(ii));
    phase1.slice(ii) =  exp(arma::cx_double(0,1)*thickness(ii)*kiz.slice(ii));

 }

  // now recursion, r.tmp is the combined reflection from last layer
  // r_{N-1,N}, then r_{N-2,N}, ... finally r_{1N}
  // starting from last rsingle

  r = rsingle.slice(Nsurface - 1);
  t =  tsingle.slice(Nsurface - 1);

  for (ii = Nsurface - 2; ii >= 0; ii--){
      t = (tsingle.slice(ii) % t % phase1.slice(ii+1) ) /
  	  (1 + rsingle.slice(ii) % r % phase2.slice(ii+1));
      r = (rsingle.slice(ii) + r % phase2.slice(ii+1) ) /
  	  (1 + rsingle.slice(ii) % r % phase2.slice(ii+1));
  }

  return List::create(
		      _["reflection"]  = r,
		      _["transmission"]  = t
 ) ;

}
