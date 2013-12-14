\name{gaussian_field}
\alias{gaussian_field}
\title{gaussian_field}
\usage{
gaussian_field(x = 1, y = 1, z = 1, wavelength = 500, alpha = 15 *
  pi/180, psi = 0, w0 = 10000, ni = 1.5, no = 1, nl = no, d = 0,
  cutoff = min(1, 3 * wavelength/(ni * pi * w0)), maxEval = 3000)
}
\arguments{
  \item{x}{position}

  \item{y}{position}

  \item{z}{position}

  \item{wavelength}{wavelength}

  \item{alpha}{beam incident angle}

  \item{theta}{beam polarisation angle}

  \item{w0}{beam waist radius}

  \item{ni}{incident medium (prism) refractive index}

  \item{no}{outer medium (substrate) refractive index}

  \item{nl}{layer refractive index}

  \item{d}{thickness of layer}

  \item{cutoff}{radial integration limit}

  \item{maxEval}{passed to adaptIntegrate}
}
\value{
data.frame Electric field (squared modulus) at the x, y, z
position
}
\description{
Electric field from the transmission of a gaussian beam at
a planar interface
}
\details{
Integration is performed over a spectrum of incident plane
waves using integrand_gb
}
\author{
Baptiste Auguie
}

