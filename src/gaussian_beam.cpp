// 
// near-field profile outside a multilayer
// with excitation from a gaussian beam
// 

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma ;
using namespace std;

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

// TODO
// include in main function, but take out the calculation of fresnel coefficients with an external
// routine, to switch between single / multiple interfaces

arma::cx_colvec evanescent_field(const cx_colvec& ei2p, const colvec& ki2p, \
   const cx_colvec& ko2p,  const cx_double ni, const cx_double no)
  {

    cx_double i = cx_double(0,1);
    cx_double ni2 = ni*ni, no2 = no*no, nratio2 = ni2/no2;
    cx_colvec eo2p(3);
    cx_double tp, ts;
   
    tp = 2.0 * no2 * ki2p(2) / (no2*ki2p(2) + ni2*ko2p(2));
    ts = 2.0 * ki2p(2) / (ki2p(2) + ko2p(2));

    eo2p(0) = nratio2 * tp * ko2p(2)/ki2p(2) * ei2p(0);
    eo2p(1) = ts * ei2p(1);
    eo2p(2) = nratio2 * tp * ei2p(2);

    return (eo2p);
  }



// [[Rcpp::export]]
arma::colvec integrand_gb(const colvec& rt, const colvec& r2, const double ki, \
			  const double psi, const double alpha, const double w0, \
			  const double ni, const double no, const cx_double nl, const double d)
  {

    double delta, rho, theta, sx, sy;
    cx_double i = cx_double(0,1), a, pw, ko, kl;
  
    cx_colvec ei1(3), ei2(3), ei2p(3), eo2p(3), Eo2(3);
    colvec ki1(3), ki2(3), ki2p(3);
    cx_colvec ko2(3), ko2p(3), kl2p(3);
    mat Ry(3, 3), Rz(3, 3), Rzi(3, 3);

    cx_double ni2 = ni*ni, no2 = no*no, nl2=nl*nl, nratio2 = ni2/no2;
    cx_double tp, ts, ap, bp;

    // change of variables from polar coordinates
    rho = rt(0); theta = rt(1);
    sx = rho * cos(theta);
    sy = rho * sin(theta);

    // work out kz component
    ki1(0) = ki*sx;
    ki1(1) = ki*sy;
    ki1(2) = ki*sqrt(1 - sx*sx - sy*sy); // only real freqs.

    // incident field polarisation and distribution
    ei1 = incident_field(psi);
    a = w0*w0 / (4*datum::pi) * exp(-(ki1(0)*ki1(0) + ki1(1)*ki1(1))*(w0*w0/4));

    // rotations of incident field

    // to frame F2
    Ry = rotation_y(alpha);
    ki2 = Ry * ki1;
    ei2 = Ry * ei1;
    
    // to frame F2p
    delta = asin(ki2(1) / sqrt(ki2(0)*ki2(0) + ki2(1)*ki2(1)));
    Rz = rotation_z(delta);
    Rzi = rotation_z(-delta);
    ki2p = Rz * ki2;
    ei2p = Rz * ei2;
  
    // transmission through interface
    ko = no / ni * ki; // outer medium
    kl = nl / ni *ki; // layer
    kl2p(0) = ki2p(0); 
    kl2p(1) = ki2p(1);
    kl2p(2) = sqrt(kl*kl - kl2p(0)*kl2p(0) - kl2p(1)*kl2p(1));
    ko2p(0) = ki2p(0); 
    ko2p(1) = ki2p(1);
    ko2p(2) = sqrt(ko*ko - ko2p(0)*ko2p(0) - ko2p(1)*ko2p(1));
    // eo2p = evanescent_field(ei2p, ki2p, ko2p, ni, no);

    // Fresnel coefficients

    // tp = 2.0 * no2 * ki2p(2) / (no2*ki2p(2) + ni2*ko2p(2));
    // ts = 2.0 * ki2p(2) / (ki2p(2) + ko2p(2));

    // p-pol
    ap = ki2p(2) / ni2;
    bp = kl2p(2) / nl2;
    arma::cx_double rp01 = (ap - bp) / (ap + bp);
    arma::cx_double tp01 = 2.0 * ap / (ap + bp);
    
    ap = kl2p(2) / nl2;
    bp = ko2p(2) / no2;
    arma::cx_double rp12 = (ap - bp) / (ap + bp);
    arma::cx_double tp12 = 2.0 * ap / (ap + bp);

    // s-pol

    ap = ki2p(2);
    bp = kl2p(2);
    arma::cx_double rs01 = (ap - bp) / (ap + bp);
    arma::cx_double ts01 = 2.0 * ap / (ap + bp);
    
    ap = kl2p(2);
    bp = ko2p(2);
    arma::cx_double rs12 = (ap - bp) / (ap + bp);
    arma::cx_double ts12 = 2.0 * ap / (ap + bp);
 
    arma::cx_double phase1 =  exp(i*d*kl2p(2));
    arma::cx_double phase2 =  exp(2.0*i*d*kl2p(2));
   
    // rp  = ( rp01 + rp12%phase2 ) / ( 1 + rp01%rp12%phase2 );
    // rs  = ( rs01 + rs12%phase2 ) / ( 1 + rs01%rs12%phase2 );
    tp  = ( tp01 * tp12 * phase1 ) / ( 1.0 + rp01 * rp12 * phase2 );
    ts  = ( ts01 * ts12 * phase1 ) / ( 1.0 + rs01 * rs12 * phase2 );
    
    eo2p(0) = nratio2 * tp * ko2p(2)/ki2p(2) * ei2p(0);
    eo2p(1) = ts * ei2p(1);
    eo2p(2) = nratio2 * tp * ei2p(2);

    // incident plane wave
    pw = exp(i*(ki2(0)*r2(0) + ki2(1)*r2(1) + sqrt(ko*ko - ki2(0)*ki2(0) +0.0*i)*r2(2)));

    Eo2 = rho * a * pw  * Rzi * eo2p;

    // join real and imaginary part in 6-vector for cubature::adaptIntegrate
    colvec Er = real(Eo2);
    colvec Ei = imag(Eo2);
    colvec res = join_cols(Er, Ei);

    return (res);
  }

RCPP_MODULE(gaussian){
  Rcpp::function( "integrand_gb", &integrand_gb,					\
	    "Integrand for the transmitted field under gaussian illumination" ) ;

}
