\name{multilayer}
\alias{multilayer}
\title{multilayer}
\usage{
  multilayer(lambda = NULL, k0 = 2 * pi/lambda,
    theta = NULL, q = sin(theta),
    epsilon = list(incident = 1.5^2, 1.33),
    thickness = c(0, 0), polarisation = c("p", "s"), d = 1,
    dout = d, ...)
}
\arguments{
  \item{lambda}{[vector] wavelength in nm}

  \item{k0}{[vector] wavevector in nm^-1}

  \item{theta}{[vector] incident angles in radians}

  \item{q}{[vector] normalised incident in-plane
  wavevector}

  \item{epsilon}{list of N+2 dielectric functions, each of
  length 1 or length(lambda)}

  \item{thickness}{vector of N+2 layer thicknesses, first
  and last are dummy}

  \item{d}{vector of distances where LFIEF are evaluated
  from each interface}

  \item{dout}{vector of distances where LFIEF are evaluated
  outside the stack}

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
\references{
  Principles of surface-enhanced Raman spectroscopy and
  related plasmonic effects

  Eric C. Le Ru and Pablo G. Etchegoin, published by
  Elsevier, Amsterdam (2009).
}

