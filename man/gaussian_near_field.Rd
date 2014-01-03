\name{gaussian_near_field}
\alias{gaussian_near_field}
\title{gaussian_near_field}
\usage{
gaussian_near_field(xyz, wavelength = 500, alpha = 15 * pi/180, psi = 0,
  w0 = 10000, epsilon = c(1.5^2, epsAg(wavelength)$epsilon, 1^2),
  thickness = c(0, 50, 0), maxEval = 3000, tol = 1e-04,
  progress = FALSE, field = FALSE)
}
\arguments{
  \item{xyz}{position}

  \item{wavelength}{wavelength}

  \item{alpha}{beam incident angle}

  \item{psi}{beam polarisation angle}

  \item{w0}{beam waist radius}

  \item{epsilon}{vector of permittivities}

  \item{thickness}{thickness corresponding to each medium}

  \item{maxEval}{passed to adaptIntegrate}

  \item{tol}{passed to adaptIntegrate}

  \item{progress}{logical: display progress bar}

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
waves
}
\author{
Baptiste Auguie
}
\seealso{
Other gaussian_beam: \code{\link{gaussian_near_field2}}
}

