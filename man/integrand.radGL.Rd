\name{integrand.radGL}
\alias{integrand.radGL}
\title{integrand.radGL}
\usage{
  integrand.radGL(d = 10, theta, lambda,
    epsilon = list(incident = 1.5^2, 1^2),
    thickness = c(0, 0))
}
\arguments{
  \item{d}{distance in nm}

  \item{theta}{angle in radians}

  \item{lambda}{wavelength in nm}

  \item{epsilon}{list of dielectric functions}

  \item{thickness}{list of layer thicknesses}
}
\description{
  Dipole decay rates near a multilayer interface
}
\details{
  Integrand of the radiative dipole decay rates near a
  multilayer interface
}
\author{
  baptiste Auguie
}

