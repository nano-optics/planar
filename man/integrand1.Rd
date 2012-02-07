\name{integrand1}
\alias{integrand1}
\title{integrand1}
\usage{
  integrand1(d = 10, u, lambda,
    epsilon = list(incident = 1.5^2, 1^2),
    thickness = c(0, 0))
}
\arguments{
  \item{d}{distance in nm}

  \item{u}{transformed normalised in-plane wavevector
  sqrt(1-q^2)}

  \item{lambda}{wavelength in nm}

  \item{epsilon}{list of dielectric functions}

  \item{thickness}{list of layer thicknesses}
}
\description{
  Dipole decay rates near a multilayer interface
}
\details{
  Integrand of the dipole decay rates near a multilayer
  interface. Transformed part I (radiative) from u=0 to 1
}
\author{
  baptiste Auguie
}

