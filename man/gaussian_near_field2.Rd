\name{gaussian_near_field2}
\alias{gaussian_near_field2}
\title{gaussian_near_field2}
\usage{
gaussian_near_field2(x = 1, y = 1, z = 1, wavelength = 632.8,
  alpha = 15 * pi/180, psi = 0, w0 = 10000, epsilon = c(1.5^2,
  epsAg(lambda)$epsilon, 1^2, 1^2), thickness = c(0, 50, 10, 0),
  cutoff = min(1, 3 * wavelength/(Re(sqrt(epsilon[1])) * pi * w0)),
  maxEval = 3000, tol = 1e-04, field = FALSE)
}
\arguments{
  \item{x}{position}

  \item{y}{position}

  \item{z}{position}

  \item{wavelength}{wavelength}

  \item{alpha}{beam incident angle}

  \item{psi}{beam polarisation angle}

  \item{w0}{beam waist radius}

  \item{epsilon}{vector of permittivities}

  \item{thickness}{thickness corresponding to each medium}

  \item{cutoff}{radial integration limit}

  \item{maxEval}{passed to adaptIntegrate}

  \item{tol}{passed to adaptIntegrate}

  \item{field}{logical: return the electric field (complex
  vector), or modulus squared}
}
\value{
data.frame electric field at the x, y, z position
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
\seealso{
Other gaussian_beam: \code{\link{gaussian_near_field}}
}

