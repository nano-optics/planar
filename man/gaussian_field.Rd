\name{gaussian_field}
\alias{gaussian_field}
\title{gaussian_field}
\usage{
  gaussian_field(x = 1, y = 1, z = 1, lambda = 500,
    alpha = pi/2 - 15 * pi/180, psi = 0, gamma = 10,
    np = 1.5, ns = 1)
}
\arguments{
  \item{x}{position}

  \item{y}{position}

  \item{z}{position}

  \item{lambda}{wavelength}

  \item{alpha}{beam incident angle}

  \item{psi}{beam polarisation angle}

  \item{gamma}{beam waist radius}

  \item{np}{incident medium (prism) refractive index}

  \item{ns}{outer medium (substrate) refractive index}
}
\value{
  data.frame Electric field (squared modulus) at the xyz
  position
}
\description{
  Electric field from the transmission of a gaussian beam
  at a planar interface
}
\details{
  Integration is performed over a spectrum of incident
  plane waves using integrand_gb
}
\author{
  Baptiste Auguie
}

