\name{dipole_direct}
\alias{dipole_direct}
\title{dipole_direct}
\usage{
  dipole_direct(d = 1, lambda,
    epsilon = list(incident = 1^2), thickness = c(0, 0),
    Nquadrature1 = 50, Nquadrature2 = 200,
    Nquadrature3 = 50, qcut = NULL, qmax = Inf)
}
\arguments{
  \item{d}{distance in nm}

  \item{lambda}{wavelength in nm}

  \item{epsilon}{list of dielectric functions}

  \item{thickness}{list of layer thicknesses}

  \item{Nquadrature1}{quadrature points in radiative
  region}

  \item{Nquadrature2}{quadrature points in SPPs region}

  \item{Nquadrature3}{quadrature points in dipole image
  region}

  \item{qcut}{transition between regions 2 and 3}

  \item{qmax}{maximum q of region 3}
}
\description{
  Dipole total decay rate near a multilayer interface
}
\details{
  direct application of the textbook formula; performs
  poorly compared to the transformed version in
  \code{dipole}
}
\author{
  baptiste Auguie
}

