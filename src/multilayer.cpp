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


Rcpp::List layer_fresnel(const arma::colvec& k0,			\
			     const arma::cx_mat& kx,			\
			     const arma::cx_mat& epsilon,		\
			     const double& thickness) {
  

  const int Nlambda = k0.n_rows;
  const int Ntheta = kx.n_cols;
  
  const arma::colvec& k02 = k0 % k0;
  const arma::cx_mat& kx2 = kx % kx;
  
  // kiz
  arma::cx_mat kiz0 = arma::ones<arma::cx_mat>(Nlambda, Ntheta), kiz1=kiz0, kiz2=kiz0;
  
  // cout << "here" << "\n";

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

Rcpp::List multilayer_field(const double k0,				
			    const double kx,			
			    const arma::cx_vec& epsilon,		
			    const arma::colvec& thickness, 
			    const arma::colvec& z, const double psi) {
  
  // psi is the polarisation angle
  
  const int Nlayer = thickness.n_elem;

  const double k02 = k0 * k0;
  const arma::cx_double q = kx / (sqrt(epsilon(0))*k0);
  const arma::cx_double kx2 = kx * kx;
  const arma::cx_double I = arma::cx_double(0,1), oI = arma::cx_double(0,0);

  // loop to calculate kiz
  arma::cx_vec kiz = arma::zeros<arma::cx_vec>(Nlayer);
  int ii=0;
  for(ii=0; ii<Nlayer; ii++){
    kiz(ii)  = sqrt(epsilon(ii) * k02 - kx2);
  } 

  arma::cx_double Mp11=1.0, Mp12=0.0, Mp21=0.0, Mp22=1.0, M11new, M12new, M21new, M22new;
  arma::cx_double Ms11=1.0, Ms12=0.0, Ms21=0.0, Ms22=1.0;
  arma::cx_vec Mpi11=arma::ones<arma::cx_vec>(Nlayer-1), Mpi12=Mpi11, Mpi21=Mpi11, Mpi22=Mpi11;
  arma::cx_vec Msi11=arma::ones<arma::cx_vec>(Nlayer-1), Msi12=Msi11, Msi21=Msi11, Msi22=Msi11;

  arma::cx_double Kpi, Ksi, phasei;

  // loop from first to last interface, transition matrices
  for(ii=0; ii<Nlayer-1; ii++){
    
    Ksi = kiz(ii+1)/kiz(ii) ;
    Kpi = Ksi * epsilon(ii)/epsilon(ii+1);
    
    phasei = exp( I * thickness(ii) * kiz(ii)) ;
    
    Mpi11(ii) = 0.5*(1.0 +Kpi) / phasei;
    Mpi21(ii) = 0.5*(1.0 -Kpi) * phasei;
    Mpi12(ii) = 0.5*(1.0 -Kpi) / phasei;
    Mpi22(ii) = 0.5*(1.0 +Kpi) * phasei;
    M11new = Mp11*Mpi11(ii) + Mp12*Mpi21(ii);
    M21new = Mp21*Mpi11(ii) + Mp22*Mpi21(ii);
    M12new = Mp11*Mpi12(ii) + Mp12*Mpi22(ii);
    M22new = Mp21*Mpi12(ii) + Mp22*Mpi22(ii);
    Mp11 = M11new;
    Mp12 = M12new;
    Mp21 = M21new;
    Mp22 = M22new;

    Msi11(ii) = 0.5*(1.0 +Ksi) / phasei;
    Msi21(ii) = 0.5*(1.0 -Ksi) * phasei;
    Msi12(ii) = 0.5*(1.0 -Ksi) / phasei;
    Msi22(ii) = 0.5*(1.0 +Ksi) * phasei;
    M11new = Ms11*Msi11(ii) + Ms12*Msi21(ii);
    M21new = Ms21*Msi11(ii) + Ms22*Msi21(ii);
    M12new = Ms11*Msi12(ii) + Ms12*Msi22(ii);
    M22new = Ms21*Msi12(ii) + Ms22*Msi22(ii);
    Ms11 = M11new;
    Ms12 = M12new;
    Ms21 = M21new;
    Ms22 = M22new;

  } 
  arma::cx_double ts = 1.0 / Ms11, tp = 1.0 / Mp11;
  arma::cx_double rs = Ms21 * ts, rp = Mp21 * tp;
   
  // calculate the fields
  arma::cx_vec HiyH1y(Nlayer) , HpiyH1y(Nlayer) , EixE1(Nlayer) , EpixE1(Nlayer) , EizE1(Nlayer) , EpizE1(Nlayer);
  arma::cx_vec EiyE1y(Nlayer) , EpiyE1y(Nlayer) , HixH1(Nlayer) , HpixH1(Nlayer) , HizH1(Nlayer) , HpizH1(Nlayer);
  arma::cx_double AuxE1, AuxE2, AuxH1, AuxH2;


  // p-pol
  HiyH1y(Nlayer-1) = tp;
  HpiyH1y(Nlayer-1) = oI;

  AuxE1 = sqrt(epsilon(0)) / k0 / epsilon(Nlayer-1);
  EixE1(Nlayer-1) = HiyH1y(Nlayer-1) * kiz(Nlayer-1) * AuxE1;
  EpixE1(Nlayer-1) = oI;
  AuxE2 = epsilon(0) / epsilon(Nlayer-1) * q.real();
  EizE1(Nlayer-1) = - HiyH1y(Nlayer-1) * AuxE2;
  EpizE1(Nlayer-1) = oI   ;

  // s-polarisation
  EiyE1y(Nlayer-1) = ts;
  EpiyE1y(Nlayer-1) = oI;

  AuxH1 = 1.0 / (k0 * sqrt(epsilon(0)));
  HixH1(Nlayer-1) = - EiyE1y(Nlayer-1) * kiz(Nlayer-1) * AuxH1;
  HpixH1(Nlayer-1) = oI;
  AuxH2 = q.real();
  HizH1(Nlayer-1) = EiyE1y(Nlayer-1) * AuxH2;
  HpizH1(Nlayer-1) = oI;
    

  // loop downwards to compute all field amplitudes
  for(ii=Nlayer-2; ii>=0; ii--){
    // p-pol
    HiyH1y(ii) = Mpi11(ii)*HiyH1y(ii+1) + Mpi12(ii)*HpiyH1y(ii+1);
    HpiyH1y(ii) = Mpi21(ii)*HiyH1y(ii+1) + Mpi22(ii)*HpiyH1y(ii+1);

    AuxE1 = sqrt(epsilon(0)) / k0 / epsilon(ii);
    EixE1(ii) = HiyH1y(ii) * kiz(ii) * AuxE1;
    EpixE1(ii) = - HpiyH1y(ii) * kiz(ii) * AuxE1;
    AuxE2 = epsilon(0) / epsilon(ii) * q.real();
    EizE1(ii) = - HiyH1y(ii) * AuxE2;
    EpizE1(ii) = - HpiyH1y(ii) * AuxE2;

    // s-pol
    EiyE1y(ii) = Msi11(ii)*EiyE1y(ii+1) + Msi12(ii)*EpiyE1y(ii+1);
    EpiyE1y(ii) = Msi21(ii)*EiyE1y(ii+1) + Msi22(ii)*EpiyE1y(ii+1);

    HixH1(ii)  = - EiyE1y(ii) * kiz(ii) * AuxH1;
    HpixH1(ii) =   EpiyE1y(ii) * kiz(ii) * AuxH1;
    HizH1(ii)  =   EiyE1y(ii) * AuxH2;
    HpizH1(ii) =   EpiyE1y(ii) * AuxH2;

  }


  // field at distances z
  arma::cx_mat Einternal(3, z.n_elem); 
  int jj;
  for (jj=0; jj<z.n_elem; jj++) {

    arma::colvec interfaces = cumsum(thickness);
    int position; // test in which layer we are
    double d, ztmp=z(jj); // current z-position (absolute and relative)
    
    if (ztmp < 0) {
      position = 0;
      d = - ztmp;
    } else if (ztmp > interfaces(Nlayer-1)){
      position = Nlayer-1;
    d = ztmp - interfaces(Nlayer-2);
    } else {
      for (ii=0; ii<Nlayer-1; ii++) {
	if(ztmp >= interfaces(ii) & ztmp < interfaces(ii+1)) position = ii+1;
	d = ztmp - interfaces(position-1);
    }
    }
    
    // note: normalise the field components to the incident field projection along s- and p-polarisations
    // p-pol corresponds to psi=0, s-pol to psi=pi/2
    
    Einternal(0,jj) = cos(psi)*(EixE1(position)  * exp(I*d*kiz(position)) + EpixE1(position)  * exp(-I*d*kiz(position)));
    Einternal(1,jj) = sin(psi)*(EiyE1y(position) * exp(I*d*kiz(position)) + EpiyE1y(position) * exp(-I*d*kiz(position)));
    Einternal(2,jj) = cos(psi)*(EizE1(position)  * exp(I*d*kiz(position)) + EpizE1(position)  * exp(-I*d*kiz(position)));
  }

   return List::create( 
		       _["rs"]  = rs, 
		       _["rp"]  = rp, 
		       _["ts"]  = ts,
		       _["tp"]  = tp,
		       _["E"]  = Einternal ) ;
   
}

Rcpp::List multilayer(const arma::colvec& k0,				\
		      const arma::cx_mat& kx,				\
		      const arma::cx_mat& epsilon,			\
		      const arma::colvec& thickness,			\
		      const int& polarisation) {
  
  const int Nlambda = k0.n_rows;
  const int Ntheta = kx.n_cols;
  const int Nlayer = thickness.n_rows;

  const arma::colvec& k02 = k0 % k0;
  const arma::cx_mat& kx2 = kx % kx;
  const arma::cx_double I = arma::cx_double(0,1);

  // loop to calculate kiz
  arma::cx_cube kiz = arma::ones<arma::cx_cube>(Nlambda, Ntheta, Nlayer);

  // loop from first to last medium, define normal kiz
  int ii=0;
  for(ii=0; ii<Nlayer; ii++){
    kiz.slice(ii)  = sqrt(repmat(epsilon.col(ii) % k02, 1, Ntheta)  - kx2);
  } 

  arma::cx_mat M11=arma::cx_mat(Nlambda, Ntheta);
  arma::cx_mat M12=M11, M21=M11, M22=M11, M11new=M11, M12new=M11, M21new=M11, M22new=M11;
  M11.fill(1), M12.fill(0), M21.fill(0), M22.fill(1);

  arma::cx_cube Mi11=arma::ones<arma::cx_cube>(Nlambda,Ntheta,Nlayer-1), Mi12=Mi11, Mi21=Mi11, Mi22=Mi11;
  arma::cx_mat Ki= arma::cx_mat(Nlambda, Ntheta), phasei = arma::cx_mat(Nlambda, Ntheta);
  
  // loop from first to last interface (== medium - 1), transition matrices
  for(ii=0; ii<Nlayer-1; ii++){
    
    if (polarisation == 0){ // p-polarisation
      Ki = repmat(epsilon.col(ii) / epsilon.col(ii+1), 1, Ntheta) %	\
	kiz.slice(ii+1) / kiz.slice( ii) ;
    } else { // s-polarisation
      Ki = kiz.slice(ii+1) / kiz.slice( ii) ;
    }

    phasei = exp(arma::cx_double(0,1)*thickness(ii)*kiz.slice(ii)) ;
    
    Mi11.slice(ii) = 0.5*(1+Ki) / phasei;
    Mi21.slice(ii) = 0.5*(1-Ki) % phasei;
    Mi12.slice(ii) = 0.5*(1-Ki) / phasei;
    Mi22.slice(ii) = 0.5*(1+Ki) % phasei;

    M11new = M11 % Mi11.slice(ii) + M12 % Mi21.slice(ii) ;
    M21new = M21 % Mi11.slice(ii) + M22 % Mi21.slice(ii) ;
    M12new = M11 % Mi12.slice(ii) + M12 % Mi22.slice(ii) ;
    M22new = M21 % Mi12.slice(ii) + M22 % Mi22.slice(ii) ;
    
   M11 = M11new;
   M12 = M12new;
   M21 = M21new;
   M22 = M22new;

  } 
  arma::cx_mat transmission = 1 / M11;
  arma::cx_mat reflection = M21 % transmission;
    
  return List::create( 
   _["reflection"]  = reflection, 
   _["transmission"]  = transmission
 ) ;

}


Rcpp::List recursive_fresnel(const arma::colvec& k0,			\
			     const arma::cx_mat& kx,			\
			     const arma::cx_mat& epsilon,		\
			     const arma::colvec& thickness,		\
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
      t = (tsingle.slice(ii) % t % phase1.slice(ii+1) ) / \
  	  (1 + rsingle.slice(ii) % r % phase2.slice(ii+1));
      r = (rsingle.slice(ii) + r % phase2.slice(ii+1) ) / \
  	  (1 + rsingle.slice(ii) % r % phase2.slice(ii+1));
  }
  
  return List::create( 
		      _["reflection"]  = r,
		      _["transmission"]  = t
 ) ;

}

RCPP_MODULE(planar){
  Rcpp::function( "layer_fresnel", &layer_fresnel,					\
	    "Calculates the reflection and transmission coefficients of a single layer" ) ;

  Rcpp::function( "multilayer_field", &multilayer_field,					\
	    "Calculates the reflection and transmission coefficients of a multilayer stack" ) ;

  Rcpp::function( "multilayer", &multilayer,					\
	    "Calculates the reflection and transmission coefficients of a multilayer stack" ) ;

  Rcpp::function( "recursive_fresnel", &recursive_fresnel,					\
	    "Calculates the reflection coefficient of a multilayer stack" ) ;

}
