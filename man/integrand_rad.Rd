\name{integrand_rad}
\alias{integrand_rad}
\title{integrand_rad}
\usage{
  integrand_rad(d = 10, theta, lambda,
    epsilon = list(incident = 1.5^2, 1^2),
    thickness = c(0, 0), GL = FALSE)
}
\arguments{
  \item{d}{distance in nm}

  \item{theta}{angle in radians}

  \item{lambda}{wavelength in nm}

  \item{epsilon}{list of dielectric functions}

  \item{thickness}{list of layer thicknesses}

  \item{GL}{logical: result formatted for use with Gauss
  Legendre quadrature}
}
\description{
  Dipole decay rates near a multilayer interface
}
\details{
  Integrand of the radiative dipole decay rates near a
  multilayer interface.
}
\author{
  baptiste Auguie
}
\seealso{
  Other integrands dipole: \code{\link{integrand_nr1}},
  \code{\link{integrand_nr2}}, \code{\link{integrand_nr3}}
}

