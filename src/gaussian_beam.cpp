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


RCPP_MODULE(gaussian){

  Rcpp::function( "integrand_gb", &integrand_gb,			\
		  "Integrand for the transmitted field under gaussian illumination" ) ;
  
}
