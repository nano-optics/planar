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

    phasei = exp(I*thickness(ii)*kiz.slice(ii)) ;
    
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
  const arma::cx_double I = arma::cx_double(0,1);
  
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
  arma::cx_mat a = arma::cx_mat(Nlambda, Ntheta), b=a;

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
    tsingle.slice(ii) = 2 * a / (a + b);
   
    phase2.slice(ii) =  exp(arma::cx_double(0,2)*thickness(ii)*kiz.slice(ii));
    phase1.slice(ii) =  exp(arma::cx_double(0,1)*thickness(ii)*kiz.slice(ii));
   
 }
  // now recursion, r.tmp is the combined reflection
  // {r12, 23, ...}, then {r13, r24, ...}
  // finally {r1N}
  // starting from rsingle

  arma::cx_cube r = rsingle;
  arma::cx_cube t = tsingle;

  int jj=0;
  for (ii = 1; ii < Nsurface; ii++){
    for (jj = 0; jj < Nlayer - ii - 1; jj++){
      // cout << "ii=" << ii << ",jj=" << jj << "\n";
      r.slice(jj)  = (rsingle.slice(jj) + r.slice(jj+1)%phase2.slice(jj+1) ) / \
  	             (1 + rsingle.slice(jj)%r.slice(jj+1)%phase2.slice(jj+1));
      t.slice(jj)  = (tsingle.slice(jj) % t.slice(jj+1)%phase1.slice(jj+1) ) / \
  	             (1 + rsingle.slice(jj)%r.slice(jj+1)%phase2.slice(jj+1));
    }
  }
  
  

  return List::create( 
		      _["reflection"]  = r.slice(0),
		      _["transmission"]  = t.slice(0)
 ) ;

}

RCPP_MODULE(planar){
  using namespace Rcpp ;
  
  function( "multilayer", &multilayer,					\
	    "Calculates the reflection and transmission coefficients of a multilayer stack" ) ;

  function( "recursive_fresnel", &recursive_fresnel,					\
	    "Calculates the reflection coefficient of a multilayer stack" ) ;

}
