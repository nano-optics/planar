\name{dipole}
\alias{dipole}
\title{dipole2}
\usage{
  dipole(d = 1, lambda, epsilon = list(incident = 1^2),
    thickness = c(0, 0), qcut = NULL, rel.err = 0.001,
    Nquadrature1 = 1000, Nquadrature2 = 10000,
    Nquadrature3 = 10000, show.messages = TRUE)
}
\arguments{
  \item{d}{distance in nm}

  \item{lambda}{wavelength in nm}

  \item{epsilon}{list of dielectric functions}

  \item{thickness}{list of layer thicknesses}

  \item{Nquadrature1}{maximum number of quadrature points
  in radiative region}

  \item{Nquadrature2}{maximum number of quadrature points
  in SPPs region}

  \item{Nquadrature3}{maximum number of quadrature points
  in dipole image region}

  \item{rel.err}{relative error}

  \item{show.messages}{logical, display progress bar}

  \item{qcut}{transition between regions 2 and 3}
}
\description{
  Dipole decay rates near a multilayer interface
}
\details{
  dipole decay rates near a multilayer interface
}
\author{
  baptiste Auguie
}

