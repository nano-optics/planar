\name{integrand_nr3}
\alias{integrand_nr3}
\title{integrand_nr3}
\usage{
  integrand_nr3(d = 10, u, ucut, lambda,
    epsilon = list(incident = 1.5^2, 1^2),
    thickness = c(0, 0), GL = FALSE)
}
\arguments{
  \item{d}{distance in nm}

  \item{u}{transformed normalised in-plane wavevector
  sqrt(q^2 - 1)}

  \item{ucut}{limit of the integral}

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
  Integrand of the dipole decay rates near a multilayer
  interface. Transformed part III from u=ucut to infinity
}
\author{
  baptiste Auguie
}
\seealso{
  Other integrands dipole: \code{\link{integrand_mtot}},
  \code{\link{integrand_nr1}}, \code{\link{integrand_nr2}},
  \code{\link{integrand_rad}}
}

