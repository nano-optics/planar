// [[Rcpp::depends(RcppArmadillo)]]

//
// collection of light emitted from a multilayer
// within finite solid angle
// we use reciprocity, ie solve the problem of excitation by multiple
// plane waves and add the field intensities (incoherent)
//

#include <RcppArmadillo.h>
#include <iostream>

#include "multilayer.h"
#include "incidence.h"
#include "utils.h"
#include "cubature.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma ;
using namespace std;



//
// Integrand for collection over a range of scattering angles
//
// rt
// r2
// k0
// psi
// epsilon
// thickness
// returns a scalar intensity
// [[Rcpp::export]]
double integrand_collection(const arma::colvec& rt,
			    const arma::colvec& r2,
			    const double k0,
			    const double psi,
			    const arma::cx_vec& epsilon,
			    const arma::vec& thickness)
  {
    const int Nlayer = epsilon.n_elem;
    double delta, rho, theta, sx, sy;
    const arma::cx_double i = arma::cx_double(0,1);
    arma::cx_double a, pw;

    arma::cx_colvec Eo2 = arma::zeros<arma::cx_colvec>(3); // result
    double res = 0.0;

    // intermediate vectors
    arma::cx_colvec ei1(3), ei2(3), ei2p(3);
    arma::colvec ki1(3), ki2(3), ki2p(3);
    arma::cx_colvec ko2(3), ko2p(3);
    // rotation matrices
    arma::mat Ry(3, 3), Rz(3, 3), Rzi(3, 3);

    // define convenient wavenumbers
    arma::cx_double ni =  sqrt(epsilon(0));
    arma::cx_double no =  sqrt(epsilon(Nlayer-1));
    double ki = real(ni)*k0; // nonabsorbing incident medium
    arma::cx_double ko = no * k0; // outer medium, can be absorbing

    // change of variables from polar coordinates
    rho = sin(rt(0)); theta = rt(1);
    sx = rho * cos(theta);
    sy = rho * sin(theta);

    double root = 1.0 - rho*rho;

    // work out kz component
    ki1(0) = ki*sx;
    ki1(1) = ki*sy;
    double kpar2 = ki1(0)*ki1(0) + ki1(1)*ki1(1);
    ki1(2) = sqrt(ki*ki - kpar2); // only real freqs.
    colvec s1(3);
    s1(0) = sx;
    s1(1) = sy;
    s1(2) = sqrt(root);
    // incident field polarisation
    ei1 = incident_field2(psi, s1);

    // rotation of incident field
    // to frame F2
    // useless, we're already in F2
    //Ry = rotation_y(alpha);
    ki2 = ki1;
    ei2 = ei1;

    // outer medium
    ko2(0) = ki2(0);
    ko2(1) = ki2(1);
    //not as above, because of rotation!
    kpar2 = ki2(0)*ki2(0) + ki2(1)*ki2(1);
    double kpar = sqrt(kpar2);
    ko2(2) = sqrt(ko*ko - kpar2);

    // to frame F2p
    delta = 0.0; // case of exact normal incidence (no k-parallel)
    if (abs(kpar) >= 2*datum::eps) {
      double sindelta = ki2(1) / kpar;// safe division
      delta = asin(sindelta);
    }
    Rz = rotation_z(delta);
    Rzi = rotation_z(-delta);
    ki2p = Rz * ki2;
    ei2p = Rz * ei2;

    // wave vector on the outer side
    ko2p(0) = ki2p(0);
    ko2p(1) = ki2p(1);
    ko2p(2) = sqrt(ko*ko - ko2p(0)*ko2p(0) - ko2p(1)*ko2p(1));

    // Fresnel coefficients
    double kx = ki2p(0); // in the 2p frame, x is in the plane of incidence
    colvec z(1); // multilayer_field expects a vector
    z(0) = r2(2);

    // use psi - delta to keep E fixed along x2p (to check)
    Rcpp::List solution = cpp_multilayer_field(k0, kx, epsilon, thickness,
					   z, psi - delta);

    // cx_double rp = solution["rp"] ;
    // cx_double rs = solution["rs"] ;
    arma::cx_colvec eo2p = solution["E"] ;
    // Rzi * eo2p is the field rotated back into the fixed R2 frame
    arma::cx_colvec eo2 =  Rzi * eo2p;

    // in-plane component of plane wave
    pw = exp(i*(ko2(0)*r2(0) + ko2(1)*r2(1)));

    // rho*ki*ki from Jacobian
    // pw is the phase factor for in-plane propagation
    // eo2 is the field rotated back into the fixed R2 frame

    //Eo2 = rho*ki*ki/s1(2) * pw * eo2;
    Eo2 = pw * eo2;

    res = real(cdot(Eo2, Eo2));
    return (res);
  }

  struct collparams {
    arma::colvec r2;
    double k0;
    double psi;
    arma::cx_colvec epsilon;
    arma::colvec thickness;
  } ;


/* wrapper of integrand for integration */
int fwrapcoll(unsigned ndim, const double *x, void *fdata,
	      unsigned fdim, double *fval) {
  collparams params = *((collparams *) fdata);

  arma::colvec res(fval, 1, false);
  arma::colvec xx(2);
  xx[0] = x[0]; xx[1]=x[1];

  res = integrand_collection(xx, params.r2, params.k0, params.psi,
			     params.epsilon, params.thickness);

  return 0;
}

//
// Integrand for collection over a range of scattering angles
//
// r2
// k0
// psi
// omega
// epsilon
// thickness
// maxEval
// reqAbsError
// tol
// progress
// returns a vector of intensities
// [[Rcpp::export]]
arma::vec cpp_field_collection(const arma::mat& r2, const double k0,
			      const double psi, const arma::vec& omega,
			      const arma::cx_vec& epsilon, const arma::vec& thickness,
			      const int maxEval, const double reqAbsError,
			      const double tol, bool progress)
  {

    const int ndim = 2;
    const int fdim = 1;
    const int N = r2.n_rows;
    // initialise the vectors to store integration results
    double integral;
    std::vector<double> error(fdim);
    double* integral_pt = &integral;
    double* error_pt = &error[0];

    double xmin[2] = {omega[0],0}, xmax[2] = {omega[1],2*datum::pi};

    collparams params;
    params.k0=k0;
    params.psi=psi;
    params.epsilon=epsilon;
    params.thickness=thickness;

    arma::vec result(N);

    // initialise an Armadillo vector to use external memory
    int ii;
    for (ii=0; ii< N; ii++){
      if(progress){
	progress_bar(ii,N);
      }
      params.r2 = arma::strans(r2.row(ii));

   /* int hcubature(unsigned fdim, integrand f, void *fdata, */
   /*               unsigned dim, const double *xmin, const double *xmax,  */
   /*               size_t maxEval, double reqAbsError, double reqRelError,  */
   /*               error_norm norm, */
   /*               double *val, double *err); */

      hcubature(fdim, fwrapcoll, &params, ndim, xmin, xmax, maxEval, reqAbsError, tol, ERROR_INDIVIDUAL, integral_pt, error_pt);
      result(ii) = integral;
    }

    if(progress)
      Rcpp::Rcout << "\n";

    return (result);
  }
