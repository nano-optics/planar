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



// @title Transfer matrix formalism for multilayer stacks
// @description Full calculation for multiple wavelengths
// @param k0 vector of wavenumbers
// @param kx complex matrix of incident k parallel
// @param epsilon complex matrix of epsilon
// @param thickness vector of thicknesses
// @param z vector of positions
// @param psi scalar polarisation angle
// @param intensity boolean: also return intensities
// @return list with rs, rp, ts, tp, I (optional)
// @describeIn multilayer full multi
// @family multilayer
// [[Rcpp::export]]
Rcpp::List cpp_multilayer(const arma::colvec& k0,
			  const arma::cx_mat& kx,
			  const arma::cx_mat& epsilon,
			  const arma::colvec& thickness,
			  const arma::colvec& z,
			  const double psi,
			  const bool intensity) {

  const int Nlambda = k0.n_rows;
  const int Ntheta = kx.n_cols;
  const int Nlayer = thickness.n_elem;

  const arma::colvec k02 = k0 % k0;
  const arma::cx_mat q = kx /repmat(sqrt(epsilon.col(0)) % k0, 1, Ntheta);
  const arma::cx_mat kx2 = kx % kx;
  const arma::cx_double I = arma::cx_double(0,1);
  const arma::cx_mat oI = 0.0*arma::zeros<arma::cx_mat>(Nlambda,Ntheta);

  // loop to calculate kiz
  // loop from first to last medium, define normal kiz
  arma::cx_cube kiz = arma::ones<arma::cx_cube>(Nlambda, Ntheta, Nlayer);
  int ii=0;
  for(ii=0; ii<Nlayer; ii++){
    kiz.slice(ii)  = sqrt(repmat(epsilon.col(ii) % k02, 1, Ntheta)  - kx2);
  }

  arma::cx_mat Mp11=arma::cx_mat(Nlambda, Ntheta), Ms11=Mp11;
  arma::cx_mat M11new=Mp11, M12new=Mp11, M21new=Mp11, M22new=Mp11;
  arma::cx_mat Mp12=Mp11, Mp21=Mp11, Mp22=Mp11, Ms12=Ms11, Ms21=Ms11, Ms22=Ms11;
  //initialise
  Mp11.fill(1.0), Mp12.fill(0.0), Mp21.fill(0.0), Mp22.fill(1.0);
  Ms11.fill(1.0), Ms12.fill(0.0), Ms21.fill(0.0), Ms22.fill(1.0);

  arma::cx_cube Mpi11=arma::ones<arma::cx_cube>(Nlambda,Ntheta,Nlayer-1), Msi11=Mpi11;
  arma::cx_cube  Mpi12=Mpi11, Mpi21=Mpi11, Mpi22=Mpi11;
  arma::cx_cube  Msi12=Msi11, Msi21=Msi11, Msi22=Msi11;
  arma::cx_mat Kpi= arma::cx_mat(Nlambda, Ntheta);
  arma::cx_mat Ksi= arma::cx_mat(Nlambda, Ntheta);
  arma::cx_mat phasei = arma::cx_mat(Nlambda, Ntheta);

  // loop from first to last interface, transition matrices
  for(ii=0; ii<Nlayer-1; ii++){

      Ksi = kiz.slice(ii+1) / kiz.slice(ii) ;
      Kpi = repmat(epsilon.col(ii) / epsilon.col(ii+1), 1, Ntheta) %	\
	kiz.slice(ii+1) / kiz.slice(ii) ;

    phasei = exp(I*thickness(ii)*kiz.slice(ii)) ;
    //p-pol
    Mpi11.slice(ii) = 0.5*(1+Kpi) / phasei;
    Mpi21.slice(ii) = 0.5*(1-Kpi) % phasei;
    Mpi12.slice(ii) = 0.5*(1-Kpi) / phasei;
    Mpi22.slice(ii) = 0.5*(1+Kpi) % phasei;

    M11new = Mp11 % Mpi11.slice(ii) + Mp12 % Mpi21.slice(ii) ;
    M21new = Mp21 % Mpi11.slice(ii) + Mp22 % Mpi21.slice(ii) ;
    M12new = Mp11 % Mpi12.slice(ii) + Mp12 % Mpi22.slice(ii) ;
    M22new = Mp21 % Mpi12.slice(ii) + Mp22 % Mpi22.slice(ii) ;

    Mp11 = M11new;
    Mp12 = M12new;
    Mp21 = M21new;
    Mp22 = M22new;

    //s-pol
    Msi11.slice(ii) = 0.5*(1+Ksi) / phasei;
    Msi21.slice(ii) = 0.5*(1-Ksi) % phasei;
    Msi12.slice(ii) = 0.5*(1-Ksi) / phasei;
    Msi22.slice(ii) = 0.5*(1+Ksi) % phasei;

    M11new = Ms11 % Msi11.slice(ii) + Ms12 % Msi21.slice(ii) ;
    M21new = Ms21 % Msi11.slice(ii) + Ms22 % Msi21.slice(ii) ;
    M12new = Ms11 % Msi12.slice(ii) + Ms12 % Msi22.slice(ii) ;
    M22new = Ms21 % Msi12.slice(ii) + Ms22 % Msi22.slice(ii) ;

    Ms11 = M11new;
    Ms12 = M12new;
    Ms21 = M21new;
    Ms22 = M22new;

  }
  arma::cx_mat ts = 1.0 / Ms11, tp = 1.0 / Mp11;
  arma::cx_mat rs = Ms21 % ts, rp = Mp21 % tp;

if(intensity){
  // calculate the fields
  arma::cx_cube HiyH1y = arma::zeros<arma::cx_cube>(Nlambda,Ntheta,Nlayer),
    HpiyH1y=HiyH1y, EixE1=HiyH1y , EpixE1=HiyH1y , EizE1=HiyH1y , EpizE1=HiyH1y;
  arma::cx_cube EiyE1y = arma::zeros<arma::cx_cube>(Nlambda,Ntheta,Nlayer) ,
    EpiyE1y=EiyE1y, HixH1=EiyE1y, HpixH1=EiyE1y, HizH1=EiyE1y, HpizH1=EiyE1y;
  arma::cx_mat AuxE1(Nlambda,Ntheta), AuxE2(Nlambda,Ntheta),
    AuxH1(Nlambda,Ntheta), AuxH2(Nlambda,Ntheta);

  // p-pol
  HiyH1y.slice(Nlayer-1) = tp;
  HpiyH1y.slice(Nlayer-1) = oI;

  AuxE1 = repmat(sqrt(epsilon.col(0)) / k0 / epsilon.col(Nlayer-1), 1, Ntheta);
  EixE1.slice(Nlayer-1) = HiyH1y.slice(Nlayer-1) % kiz.slice(Nlayer-1) % AuxE1;
  EpixE1.slice(Nlayer-1) = oI;
  AuxE2 = repmat(epsilon.col(0) / epsilon.col(Nlayer-1), 1, Ntheta) % real(q);
  EizE1.slice(Nlayer-1) = - HiyH1y.slice(Nlayer-1) % AuxE2;
  EpizE1.slice(Nlayer-1) = oI   ;

  // s-polarisation
  EiyE1y.slice(Nlayer-1) = ts;
  EpiyE1y.slice(Nlayer-1) = oI;

  AuxH1 = repmat(1.0 / (sqrt(epsilon.col(0)) % k0), 1, Ntheta);
  HixH1.slice(Nlayer-1) = - EiyE1y.slice(Nlayer-1) % kiz.slice(Nlayer-1) % AuxH1;
  HpixH1.slice(Nlayer-1) = oI;
  AuxH2 = oI +real(q);
  HizH1.slice(Nlayer-1) = EiyE1y.slice(Nlayer-1) % AuxH2;
  HpizH1.slice(Nlayer-1) = oI;

  // loop downwards to compute all field amplitudes
  for(ii=Nlayer-2; ii>=0; ii--){
    // p-pol
    HiyH1y.slice(ii) = Mpi11.slice(ii) % HiyH1y.slice(ii+1) +
      Mpi12.slice(ii) % HpiyH1y.slice(ii+1);
    HpiyH1y.slice(ii) = Mpi21.slice(ii) % HiyH1y.slice(ii+1) +
      Mpi22.slice(ii) % HpiyH1y.slice(ii+1);

    AuxE1 = repmat(sqrt(epsilon.col(0)) / k0 / epsilon.col(ii), 1, Ntheta);
    EixE1.slice(ii) = HiyH1y.slice(ii) % kiz.slice(ii) % AuxE1;
    EpixE1.slice(ii) = - HpiyH1y.slice(ii) % kiz.slice(ii) % AuxE1;
    AuxE2 = repmat(epsilon.col(0) / epsilon.col(ii), 1, Ntheta) % real(q);
    EizE1.slice(ii) = - HiyH1y.slice(ii) % AuxE2;
    EpizE1.slice(ii) = - HpiyH1y.slice(ii) % AuxE2;

    // s-pol
    EiyE1y.slice(ii) = Msi11.slice(ii) % EiyE1y.slice(ii+1) +
      Msi12.slice(ii) % EpiyE1y.slice(ii+1);
    EpiyE1y.slice(ii) = Msi21.slice(ii) % EiyE1y.slice(ii+1) +
      Msi22.slice(ii) % EpiyE1y.slice(ii+1);

    HixH1.slice(ii)  = - EiyE1y.slice(ii) % kiz.slice(ii) % AuxH1;
    HpixH1.slice(ii) =   EpiyE1y.slice(ii) % kiz.slice(ii) % AuxH1;
    HizH1.slice(ii)  =   EiyE1y.slice(ii) % AuxH2;
    HpizH1.slice(ii) =   EpiyE1y.slice(ii) % AuxH2;

  }

  // tmp field at distance z
  arma::cx_mat Ex(Nlambda, Ntheta), Ey=Ex, Ez=Ex;
  // stored intensity
  arma::cube E2(Nlambda, Ntheta, z.n_elem);
  arma::colvec interfaces = cumsum(thickness);
  int position=0; // test in which layer we are

  int jj;
  for (jj=0; jj<z.n_elem; jj++) {

    double d, ztmp=z(jj); // current z-position (absolute and relative)
    if (ztmp <= 0.0) {
      position = 0;
      d = ztmp; // already negative btw...
    } else if (ztmp >= interfaces(Nlayer-1)){
      position = Nlayer-1;
      d = ztmp - interfaces(Nlayer-2);
    } else {
      for (ii=0; ii<Nlayer-2; ii++) {
	// note: could use sort_index_stable or find
	if((ztmp >= interfaces(ii)) && (ztmp < interfaces(ii+1))) {
	  position = ii+1;
	  d = ztmp - interfaces(position-1);
	  break;
	}
    }
    }
    // note: normalise the field components to the incident field projection
    //along s- and p-polarisations
    // p-pol corresponds to psi=0, s-pol to psi=pi/2
    Ex = cos(psi)*(EixE1.slice(position) % exp(I*d*kiz.slice(position)) +
		   EpixE1.slice(position) % exp(-I*d*kiz.slice(position)));
    Ey = sin(psi)*(EiyE1y.slice(position) % exp(I*d*kiz.slice(position)) +
		   EpiyE1y.slice(position) % exp(-I*d*kiz.slice(position)));
    Ez = cos(psi)*(EizE1.slice(position)  % exp(I*d*kiz.slice(position)) +
		   EpizE1.slice(position)  % exp(-I*d*kiz.slice(position)));
    E2.slice(jj) = abs(Ex) % abs(Ex) + abs(Ey) % abs(Ey) + abs(Ez) % abs(Ez);
  }

   return List::create(
		       _["rs"]  = rs,
		       _["rp"]  = rp,
		       _["ts"]  = ts,
		       _["tp"]  = tp,
		       _["I"] = E2) ;

} else { // intensity not requested

	   return List::create(
			       _["rs"]  = rs,
			       _["rp"]  = rp,
			       _["ts"]  = ts,
			       _["tp"]  = tp) ;
}
}


// @title Transfer matrix formalism for multilayer stacks
// @description Field calculation for one wavelength
// @return list with rs, rp, ts, tp, I
// @describeIn multilayer full multi
// @family multilayer
// [[Rcpp::export]]
Rcpp::List cpp_multilayer_field(const double k0,
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
  arma::cx_mat Hinternal(3, z.n_elem);
  arma::colvec interfaces = cumsum(thickness);
  int position=0; // test in which layer we are

  int jj;
  for (jj=0; jj<z.n_elem; jj++) {

    double d, ztmp=z(jj); // current z-position (absolute and relative)
    if (ztmp <= 0.0) {
      position = 0;
      d = ztmp; // already negative btw...
    } else if (ztmp >= interfaces(Nlayer-1)){
      position = Nlayer-1;
      d = ztmp - interfaces(Nlayer-2);
    } else {
      for (ii=0; ii<Nlayer-2; ii++) {
	// note: could use sort_index_stable or find
	if((ztmp >= interfaces(ii)) && (ztmp < interfaces(ii+1))) {
	  position = ii+1;
	  d = ztmp - interfaces(position-1);
	  break;
	}
    }
    }
    // note: normalise the field components to the incident field projection along s- and p-polarisations
    // p-pol corresponds to psi=0, s-pol to psi=pi/2
    Einternal(0,jj) = cos(psi)*(EixE1(position)  * exp(I*d*kiz(position)) + EpixE1(position)  * exp(-I*d*kiz(position)));
    Einternal(1,jj) = sin(psi)*(EiyE1y(position) * exp(I*d*kiz(position)) + EpiyE1y(position) * exp(-I*d*kiz(position)));
    Einternal(2,jj) = cos(psi)*(EizE1(position)  * exp(I*d*kiz(position)) + EpizE1(position)  * exp(-I*d*kiz(position)));

    // note: normalise the field components to the incident field projection along s- and p-polarisations
    // p-pol corresponds to psi=0, s-pol to psi=pi/2
		// HiyH1y HpiyH1y HixH1 HizH1 HpizH1 HpixH1
    Hinternal(0,jj) = cos(psi)*(HixH1(position)  * exp(I*d*kiz(position)) + HpixH1(position)  * exp(-I*d*kiz(position)));
    Hinternal(1,jj) = sin(psi)*(HiyH1y(position) * exp(I*d*kiz(position)) + HpiyH1y(position) * exp(-I*d*kiz(position)));
    Hinternal(2,jj) = cos(psi)*(HizH1(position)  * exp(I*d*kiz(position)) + HpizH1(position)  * exp(-I*d*kiz(position)));
  }

   return List::create(
		       _["rs"]  = rs,
		       _["rp"]  = rp,
		       _["ts"]  = ts,
		       _["tp"]  = tp,
		       _["E"]  = Einternal ,
		       _["H"]  = Hinternal ) ;

}
