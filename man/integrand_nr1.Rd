\name{integrand_nr1}
\alias{integrand_nr1}
\title{integrand_nr1}
\usage{
  integrand_nr1(d = 10, u, wavelength,
    epsilon = list(incident = 1.5^2, 1^2),
    thickness = c(0, 0), GL = FALSE)
}
\arguments{
  \item{d}{distance in nm}

  \item{u}{transformed normalised in-plane wavevector
  sqrt(1-q^2)}

  \item{wavelength}{wavelength in nm}

  \item{epsilon}{list of dielectric functions}

  \item{thickness}{list of layer thicknesses}

  \item{GL}{logical: result formatted for use with Gauss
  Legendre quadrature}
}
\description{
  Dipole decay rates near a multilayer interface
}
\details{
  Integrand of the dipole decay rates near a multilayer
  interface. Transformed part I1 (radiative) from u=0 to 1
}
\author{
  baptiste Auguie
}
\seealso{
  Other integrands dipole: \code{\link{integrand_mtot}},
  \code{\link{integrand_nr2}}, \code{\link{integrand_nr3}},
  \code{\link{integrand_rad}}
}

