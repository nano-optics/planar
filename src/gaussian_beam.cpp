// 
// near-field profile outside a multilayer
// with excitation from a gaussian beam
// 

#include <RcppArmadillo.h>
#include <iostream>

#include "multilayer.h"
#include "cubature.h"

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


// [[Rcpp::export]]
arma::colvec integrand_gb2(const colvec& rt, const colvec& r2, const double k0, 
			  const double psi, const double alpha, const double w0, 
			  const cx_vec& epsilon, const vec& thickness)
  {
    const int Nlayer = epsilon.n_elem;
    double delta, rho, theta, sx, sy;
    const cx_double i = cx_double(0,1);
    cx_double a, pw;

    cx_colvec Eo2 = arma::zeros<arma::cx_colvec>(3); // result
    colvec res = arma::zeros<arma::colvec>(6); // combining real and imag parts

    // intermediate vectors
    cx_colvec ei1(3), ei2(3), ei2p(3);
    colvec ki1(3), ki2(3), ki2p(3);
    cx_colvec ko2(3), ko2p(3);
    // rotation matrices
    mat Ry(3, 3), Rz(3, 3), Rzi(3, 3);

    // define convenient wavenumbers
    cx_double ni =  sqrt(epsilon(0));
    cx_double no =  sqrt(epsilon(Nlayer-1));
    double ki = real(ni)*k0; // nonabsorbing incident medium
    cx_double ko = no * k0; // outer medium, can be absorbing
    cx_double nini = epsilon(0), nono = epsilon(Nlayer-1);

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
    colvec s1(3);
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
    colvec z(1); // multilayer_field expects a vector
    z(0) = r2(2);

    // use psi - delta to keep E fixed along x2p (to check)
    Rcpp::List solution = multilayer_field(k0, kx, epsilon, thickness, z, psi - delta);
    
    // cx_double rp = solution["rp"] ;
    // cx_double rs = solution["rs"] ;
    cx_colvec eo2p = solution["E"] ;
    // Rzi * eo2p is the field rotated back into the fixed R2 frame
    cx_colvec eo2 =  Rzi * eo2p;

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

    colvec Er = real(Eo2);
    colvec Ei = imag(Eo2);
    res = join_cols(Er, Ei);

    return (res);
  }


// [[Rcpp::export]]
arma::colvec integrand_gb(const colvec& rt, const colvec& r2, const double ki, \
			  const double psi, const double alpha, const double w0, \
			  const double ni, const double no, const cx_double nl, const double d)
  {

    double delta, rho, theta, sx, sy;
    cx_double i = cx_double(0,1), a, pw, pwr, ko, kl;

    bool reflected = r2(2) < 0.0; // reflected side of the interface
  
    cx_colvec Eo2 = arma::zeros<arma::cx_colvec>(3); // result
    cx_colvec ei1(3), ei2(3), ei2p(3), eo2p(3);
    colvec ki1(3), ki2(3), ki2p(3);
    cx_colvec ko2(3), ko2p(3), kl2p(3);
    mat Ry(3, 3), Rz(3, 3), Rzi(3, 3);

    cx_double nini = ni*ni, nono = no*no, nlnl=nl*nl;
    ko = no / ni * ki; // outer medium
    kl = nl / ni *ki; // layer
    cx_double tp, ts, rp, rs, ap, bp, as, bs;
    // temporary variables for internal field calculations
    cx_double Kp, Ks, Mp11, Ms11, Mp21, Ms21, Hy, Hyp;

    // change of variables from polar coordinates
    rho = rt(0); theta = rt(1);
    sx = rho * cos(theta);
    sy = rho * sin(theta);

    double root = 1.0 - rho*rho;

    if( root < 0.0) // only prop. waves are considered
      {    
	colvec Er = real(Eo2);
	colvec Ei = imag(Eo2);
	colvec res = join_cols(Er, Ei);
	return(res); 
      }

    // work out kz component
    ki1(0) = ki*sx;
    ki1(1) = ki*sy;
    ki1(2) = sqrt(ki*ki - ki1(0)*ki1(0) - ki1(1)*ki1(1)); // only real freqs.

    // incident field polarisation and distribution
    ei1 = incident_field(psi);
    a = w0*w0 / (4*datum::pi) * exp(-(ki1(0)*ki1(0) + ki1(1)*ki1(1))*(w0*w0/4));

    // rotations of incident field

    // to frame F2
    Ry = rotation_y(alpha);
    ki2 = Ry * ki1;
    ei2 = Ry * ei1;
    ko2(0) = ki2(0);
    ko2(1) = ki2(1);
    ko2(2) = sqrt(ko*ko - ko2(0)*ko2(0) - ko2(1)*ko2(1));
    
    // to frame F2p
    delta = asin(ki2(1) / sqrt(ki2(0)*ki2(0) + ki2(1)*ki2(1)));
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

    // internal field in layer not implemented
    // TODO: factor out the field calculation using transfer matrices
    // to deal with arbitrary number of layers etc.

    // reflection from layer 
    // note: field is calculated outside of the layer, not inside
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
      // transmitted field
      eo2p(0) = nini / nono * tp * ko2p(2)/ki2p(2) * ei2p(0);
      eo2p(1) = ts * ei2p(1);
      eo2p(2) = nini / nono * tp * ei2p(2);
      
      // incident plane wave
      pw = exp(i*(ko2(0)*r2(0) + ko2(1)*r2(1) + ko2(2)*r2(2)));
      Eo2 = rho*ki*ki * a * pw  * Rzi * eo2p;
    }

    // join real and imaginary part in 6-vector 
    // for cubature::adaptIntegrate

    colvec Er = real(Eo2);
    colvec Ei = imag(Eo2);
    colvec res = join_cols(Er, Ei);

    return (res);
  }


  struct params {
    colvec r2;
    double k0;
    double psi;
    double alpha;
    double w0;
    cx_colvec epsilon;
    colvec thickness;
  } ;


/* wrapper of integrand for integration */
int fwrap(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
  params mydata = *((params *) fdata);
  
  colvec res(fval, 6, false);
  colvec xx(2);
  xx[0] = x[0]; xx[1]=x[1];

  res = integrand_gb2(xx, mydata.r2, mydata.k0, mydata.psi, mydata.alpha, mydata.w0, mydata.epsilon, mydata.thickness);
  
  return 0;
}

// [[Rcpp::export]]
arma::cx_mat gb_field(const mat& r2, const double k0, 
		      const double psi, const double alpha, const double w0, 
		      const cx_vec& epsilon, const vec& thickness, 
		      const int maxEval, const double tol, bool progress)
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

    params mydata;
    mydata.k0=k0;
    mydata.psi=psi;
    mydata.alpha=alpha;
    mydata.w0=w0;
    mydata.epsilon=epsilon;
    mydata.thickness=thickness;

    cx_mat result(3,N);

    // initialise an Armadillo vector to use external memory
    vec tmp(integral_pt, fdim, false);
    int ii;
    for (ii=0; ii< N; ii++){
      if(progress){
	progress_bar(ii,N);
      }
      mydata.r2 = strans(r2.row(ii));

   /* int hcubature(unsigned fdim, integrand f, void *fdata, */
   /*               unsigned dim, const double *xmin, const double *xmax,  */
   /*               size_t maxEval, double reqAbsError, double reqRelError,  */
   /*               error_norm norm, */
   /*               double *val, double *err); */

      hcubature(fdim, fwrap, &mydata, ndim, xmin, xmax, maxEval, 0, tol, ERROR_INDIVIDUAL, integral_pt, error_pt);
      result(0,ii) = cx_double(tmp(0), tmp(3));
      result(1,ii) = cx_double(tmp(1), tmp(4));
      result(2,ii) = cx_double(tmp(2), tmp(5));
    }

    if(progress)
      Rcpp::Rcout << "\n";

    return (result);
  }


RCPP_MODULE(gaussian){

  Rcpp::function( "integrand_gb", &integrand_gb,			\
		  "Integrand for the transmitted field under gaussian illumination" ) ;

  Rcpp::function( "integrand_gb2", &integrand_gb2,			\
		  "Integrand for the transmitted field under gaussian illumination" ) ;

  Rcpp::function( "gb_field", &gb_field,			\
		  "Near field under gaussian illumination" ) ;
  
}
