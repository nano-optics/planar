\name{multilayercpp}
\alias{multilayercpp}
\title{multilayercpp}
\usage{
multilayercpp(wavelength = 2 * pi/k0, k0 = 2 * pi/wavelength,
  angle = asin(q), q = sin(angle), epsilon = list(incident = 1.5^2, 1.33),
  thickness = c(0, 0), polarisation = c("p", "s"), ...)
}
\arguments{
  \item{wavelength}{[vector] wavelength in nm}

  \item{k0}{[vector] wavevector in nm^-1}

  \item{angle}{[vector] incident angles in radians}

  \item{q}{[vector] normalised incident in-plane
  wavevector}

  \item{epsilon}{list of N+2 dielectric functions, each of
  length 1 or length(wavelength)}

  \item{thickness}{vector of N+2 layer thicknesses, first
  and last are dummy}

  \item{polarisation}{[character] switch between p- and s-
  polarisation}

  \item{...}{unused}
}
\value{
fresnel coefficients and field profiles
}
\description{
Multilayer Fresnel coefficients
}
\details{
solves the EM problem of a multilayered interface
}
\examples{
library(planar)
demo(package="planar")
}
\author{
baptiste Auguie
}

