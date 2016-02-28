// [[Rcpp::depends(RcppArmadillo)]]

//
// near-field profile outside a multilayer
// with excitation from a gaussian beam
//

#include <RcppArmadillo.h>
#include <iostream>

#include "multilayer.h"
#include "cubature.h"
#include "incidence.h"
#include "utils.h"

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
arma::colvec cpp_integrand_gb_ml(const arma::colvec& rt, const arma::colvec& r2, const double k0,
				const double psi, const double alpha, const double w0,
				const arma::cx_vec& epsilon, const arma::vec& thickness)
  {
    const int Nlayer = epsilon.n_elem;
    double delta, rho, theta, sx, sy;
    const arma::cx_double i = arma::cx_double(0,1);
    arma::cx_double a, pw;

    arma::cx_colvec Eo2 = arma::zeros<arma::cx_colvec>(3); // result
    arma::colvec res = arma::zeros<arma::colvec>(6); // combining real and imag parts

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
    rho = rt(0); theta = rt(1);
    sx = rho * cos(theta);
    sy = rho * sin(theta);

    // only prop. waves are considered
    double root = 1.0 - rho*rho;
    if( root < 0.0)  return(res);

    // work out kz component
    ki1(0) = ki*sx;
    ki1(1) = ki*sy;
    double kpar2 = ki1(0)*ki1(0) + ki1(1)*ki1(1);
    ki1(2) = sqrt(ki*ki - kpar2); // only real freqs.
    arma::colvec s1(3);
    s1(0) = sx;
    s1(1) = sy;
    s1(2) = sqrt(root);
    // incident field polarisation
    ei1 = incident_field(psi);
    // use this expression for a focused beam only?
    // ei1 = incident_field2(psi, s1);
    // weight factor
    a = w0*w0 / (4*datum::pi) * exp(-kpar2*(w0*w0/4));

    // rotation of incident field
    // to frame F2
    Ry = rotation_y(alpha);
    ki2 = Ry * ki1;
    ei2 = Ry * ei1;

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
    arma::colvec z(1); // multilayer_field expects a vector
    z(0) = r2(2);

    // use psi - delta to keep E fixed along x2p (to check)
    Rcpp::List solution = cpp_multilayer_field(k0, kx, epsilon, thickness, z, psi - delta);

    // cx_double rp = solution["rp"] ;
    // cx_double rs = solution["rs"] ;
    arma::cx_colvec eo2p = solution["E"] ;
    // Rzi * eo2p is the field rotated back into the fixed R2 frame
    arma::cx_colvec eo2 =  Rzi * eo2p;

    // in-plane component of plane wave
    pw = exp(i*(ko2(0)*r2(0) + ko2(1)*r2(1)));

    // rho*ki*ki from Jacobian
    // a is the weight factor
    // pw is the phase factor for in-plane propagation
    // eo2 is the field rotated back into the fixed R2 frame
    // for a focused beam(?), also include 1/sz
    Eo2 = rho*ki*ki * a * pw * eo2;

    // join real and imaginary part in 6-vector
    // for cubature

    arma::colvec Er = real(Eo2);
    arma::colvec Ei = imag(Eo2);
    // res = join_cols(Er, Ei);
    // now interlace Re and Im to
    // use paired convergence test in cubature
    res[0] = Er[0]; res[1] = Ei[0];
    res[2] = Er[1]; res[3] = Ei[1];
    res[4] = Er[2]; res[5] = Ei[2];
    return (res);
  }


// [[Rcpp::export]]
arma::colvec cpp_integrand_gb_layer(const arma::colvec& rt, const arma::colvec& r2, const double ki,
			  const double psi, const double alpha, const double w0,
			  const double ni, const double no, const arma::cx_double nl, const double d)
  {

    double delta, rho, theta, sx, sy;
    arma::cx_double i = arma::cx_double(0,1), a, pw, pwr, ko, kl;
    bool reflected = r2(2) < 0.0; // reflected side of the interface

    arma::cx_colvec Eo2 = arma::zeros<arma::cx_colvec>(3); // result
    arma::colvec res = arma::zeros<arma::colvec>(6); // combining real and imag parts

    arma::cx_colvec ei1(3), ei2(3), ei2p(3), eo2p(3);
    arma::colvec ki1(3), ki2(3), ki2p(3);
    arma::cx_colvec ko2(3), ko2p(3), kl2p(3);
    arma::mat Ry(3, 3), Rz(3, 3), Rzi(3, 3);

    arma::cx_double nini = ni*ni, nono = no*no, nlnl=nl*nl;
    ko = no / ni * ki; // outer medium
    kl = nl / ni *ki; // layer
    arma::cx_double tp, ts, rp, rs, ap, bp, as, bs;

    // change of variables from polar coordinates
    rho = rt(0); theta = rt(1);
    sx = rho * cos(theta);
    sy = rho * sin(theta);

    // only prop. waves are considered
    double root = 1.0 - rho*rho;
    if( root < 0.0)  return(res);

    // work out kz component
    ki1(0) = ki*sx;
    ki1(1) = ki*sy;
    double kpar2 = ki1(0)*ki1(0) + ki1(1)*ki1(1);
    ki1(2) = sqrt(ki*ki - kpar2); // only real freqs.
    arma::colvec s1(3);
    s1(0) = sx;
    s1(1) = sy;
    s1(2) = sqrt(root);
    // incident field polarisation
    ei1 = incident_field(psi);
    // use this expression for a focused beam only?
    // ei1 = incident_field2(psi, s1);
    // weight factor
    a = w0*w0 / (4*datum::pi) * exp(-kpar2*(w0*w0/4));

   // rotation of incident field
    // to frame F2
    Ry = rotation_y(alpha);
    ki2 = Ry * ki1;
    ei2 = Ry * ei1;

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

    // transmission through interface
    kl2p(0) = ki2p(0);
    kl2p(1) = ki2p(1);
    kl2p(2) = sqrt(kl*kl - kl2p(0)*kl2p(0) - kl2p(1)*kl2p(1));
    ko2p(0) = ki2p(0);
    ko2p(1) = ki2p(1);
    ko2p(2) = sqrt(ko*ko - ko2p(0)*ko2p(0) - ko2p(1)*ko2p(1));

    // Fresnel coefficients

    // p-pol
    ap = ki2p(2) / nini;
    bp = kl2p(2) / nlnl;
    arma::cx_double rp01 = (ap - bp) / (ap + bp);
    arma::cx_double tp01 = 2.0 * ap / (ap + bp);

    ap = kl2p(2) / nlnl;
    bp = ko2p(2) / nono;
    arma::cx_double rp12 = (ap - bp) / (ap + bp);
    arma::cx_double tp12 = 2.0 * ap / (ap + bp);

    // s-pol

    as = ki2p(2);
    bs = kl2p(2);
    arma::cx_double rs01 = (as - bs) / (as + bs);
    arma::cx_double ts01 = 2.0 * as / (as + bs);

    as = kl2p(2);
    bs = ko2p(2);
    arma::cx_double rs12 = (as - bs) / (as + bs);
    arma::cx_double ts12 = 2.0 * as / (as + bs);

    arma::cx_double phase1 =  exp(i*d*kl2p(2));
    arma::cx_double phase2 =  exp(2.0*i*d*kl2p(2));

    // transmission through layer
    tp  = ( tp01 * tp12 * phase1 ) / ( 1.0 + rp01 * rp12 * phase2 );
    ts  = ( ts01 * ts12 * phase1 ) / ( 1.0 + rs01 * rs12 * phase2 );

    // Note: internal field in layer not implemented
    // this means that 0 <= r2(2) < d is meaningless
    // use the other function (transfer matrix) for full internal field
    // with arbitrary number of layers (slower)

    // reflection from layer (z < 0)
    if(reflected) {
      rp  = ( rp01 + rp12 * phase2 ) / ( 1.0 + rp01 * rp12 * phase2 );
      rs  = ( rs01 + rs12 * phase2 ) / ( 1.0 + rs01 * rs12 * phase2 );

      eo2p(0) = rp*ei2p(0);
      eo2p(1) = rs*ei2p(1);
      eo2p(2) = rp*ei2p(2);

      pw = exp(i*(ki2(0)*r2(0) + ki2(1)*r2(1) + ki2(2)*r2(2))); // incident plane wave
      pwr = exp(i*(ki2(0)*r2(0) + ki2(1)*r2(1) - ki2(2)*r2(2))); // reflected plane wave
      Eo2 = rho*ki*ki * a * (pwr*Rzi*eo2p + pw*Rzi*ei2p);
    } else {
      // transmitted field (z >= d)
      eo2p(0) = nini / nono * tp * ko2p(2)/ki2p(2) * ei2p(0);
      eo2p(1) = ts * ei2p(1);
      eo2p(2) = nini / nono * tp * ei2p(2);

      // in this single-layer configuration we are using the interface as the origin
      double position = r2(2) - d;
      // incident plane wave
      pw = exp(i*(ko2(0)*r2(0) + ko2(1)*r2(1) + ko2(2)*position));
      Eo2 = rho*ki*ki * a * pw  * Rzi * eo2p;
    }

    // join real and imaginary part in 6-vector
    // for cubature

    arma::colvec Er = real(Eo2);
    arma::colvec Ei = imag(Eo2);

    res[0] = Er[0]; res[1] = Ei[0];
    res[2] = Er[1]; res[3] = Ei[1];
    res[4] = Er[2]; res[5] = Ei[2];

    return (res);
  }


  struct parameters {
    arma::colvec r2;
    double k0;
    double psi;
    double alpha;
    double w0;
    arma::cx_colvec epsilon;
    arma::colvec thickness;
  } ;


/* wrapper of integrand for integration */
int fwrap(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
  parameters params = *((parameters *) fdata);

  arma::colvec res(fval, 6, false);
  arma::colvec xx(2);
  xx[0] = x[0]; xx[1]=x[1];

// arma::colvec integrand_gb(const colvec& rt, const colvec& r2, const double ki, \
// 			  const double psi, const double alpha, const double w0, \
// 			  const double ni, const double no, const cx_double nl, const double d)

  double ni = sqrt(real(params.epsilon[0]));
  arma::cx_double nl = sqrt(params.epsilon[1]);
  double no = sqrt(real(params.epsilon[2]));
  double d = params.thickness[1];
  double ki = params.k0 * ni;

  res = cpp_integrand_gb_layer(xx, params.r2, ki, params.psi, params.alpha, params.w0, ni, no, nl, d);

  return 0;
}

/* wrapper of integrand for integration */
int fwrap2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
  parameters params = *((parameters *) fdata);

  arma::colvec res(fval, 6, false);
  arma::colvec xx(2);
  xx[0] = x[0]; xx[1]=x[1];

  res = cpp_integrand_gb_ml(xx, params.r2, params.k0, params.psi, params.alpha, params.w0, params.epsilon, params.thickness);

  return 0;
}

// [[Rcpp::export]]
arma::cx_mat cpp_field_gb_layer(const arma::mat& r2, const double k0,
		      const double psi, const double alpha, const double w0,
		      const arma::cx_vec& epsilon, const arma::vec& thickness,
		      const int maxEval, const double reqAbsError, const double tol, bool progress)
  {

    const int ndim = 2;
    const int fdim = 6;
    const int N = r2.n_rows;
    // initialise the vectors to store integration results
    std::vector<double> integral(fdim);
    std::vector<double> error(fdim);
    double* integral_pt = &integral[0];
    double* error_pt = &error[0];

    // 3*wavelength/(ni*pi*w0)
    double cutoff = 6.0 / (k0 * sqrt(real(epsilon[0])) * w0);
    double xmin[2] = {0,0}, xmax[2] = {cutoff,2*datum::pi};

    parameters params;
    params.k0=k0;
    params.psi=psi;
    params.alpha=alpha;
    params.w0=w0;
    params.epsilon=epsilon;
    params.thickness=thickness;

    arma::cx_mat result(3,N);

    // initialise an Armadillo vector to use external memory
    arma::vec tmp(integral_pt, fdim, false);
    int ii;
    for (ii=0; ii< N; ii++){
      if(progress){
	progress_bar(ii,N);
      }
      params.r2 = strans(r2.row(ii));

   /* int hcubature(unsigned fdim, integrand f, void *fdata, */
   /*               unsigned dim, const double *xmin, const double *xmax,  */
   /*               size_t maxEval, double reqAbsError, double reqRelError,  */
   /*               error_norm norm, */
   /*               double *val, double *err); */

      hcubature(fdim, fwrap, &params, ndim, xmin, xmax, maxEval, reqAbsError, tol, ERROR_PAIRED, integral_pt, error_pt);
      result(0,ii) = cx_double(tmp(0), tmp(1));
      result(1,ii) = cx_double(tmp(2), tmp(3));
      result(2,ii) = cx_double(tmp(4), tmp(5));
    }

    if(progress)
      Rcpp::Rcout << "\n";

    return (result);
  }


// [[Rcpp::export]]
arma::cx_mat cpp_field_gb_ml(const arma::mat& r2, const double k0,
		      const double psi, const double alpha, const double w0,
		      const arma::cx_vec& epsilon, const arma::vec& thickness,
		      const int maxEval, const double reqAbsError, const double tol, bool progress)
  {

    const int ndim = 2;
    const int fdim = 6;
    const int N = r2.n_rows;
    // initialise the vectors to store integration results
    std::vector<double> integral(fdim);
    std::vector<double> error(fdim);
    double* integral_pt = &integral[0];
    double* error_pt = &error[0];

    // 3*wavelength/(ni*pi*w0)
    double cutoff = 6.0 / (k0 * sqrt(real(epsilon[0])) * w0);
    double xmin[2] = {0,0}, xmax[2] = {cutoff,2*datum::pi};

    parameters params;
    params.k0=k0;
    params.psi=psi;
    params.alpha=alpha;
    params.w0=w0;
    params.epsilon=epsilon;
    params.thickness=thickness;

    arma::cx_mat result(3,N);

    // initialise an Armadillo vector to use external memory
    arma::vec tmp(integral_pt, fdim, false);
    int ii;
    for (ii=0; ii< N; ii++){
      if(progress){
	progress_bar(ii,N);
      }
      params.r2 = strans(r2.row(ii));

   /* int hcubature(unsigned fdim, integrand f, void *fdata, */
   /*               unsigned dim, const double *xmin, const double *xmax,  */
   /*               size_t maxEval, double reqAbsError, double reqRelError,  */
   /*               error_norm norm, */
   /*               double *val, double *err); */

      hcubature(fdim, fwrap2, &params, ndim, xmin, xmax, maxEval, reqAbsError, tol, ERROR_PAIRED, integral_pt, error_pt);
      result(0,ii) = arma::cx_double(tmp(0), tmp(1));
      result(1,ii) = arma::cx_double(tmp(2), tmp(3));
      result(2,ii) = arma::cx_double(tmp(4), tmp(5));
    }

    if(progress)
      Rcpp::Rcout << "\n";

    return (result);
  }
