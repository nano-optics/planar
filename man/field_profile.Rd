\name{field_profile}
\alias{field_profile}
\title{field_profile}
\usage{
  field_profile(lambda = 500, theta = 0,
    polarisation = "p", thickness = c(0, 20, 140, 20, 0),
    dmax = 200, res = 1000,
    epsilon = list(1^2, -12, 1.38^2, -12, 1.46^2), ...)
}
\arguments{
  \item{lambda}{wavelength}

  \item{theta}{angle}

  \item{dmax}{maximum distance to interface, if > layer
  thickness}

  \item{thickness}{vector of layer thickness}

  \item{res}{resolution of sampling points}

  \item{epsilon}{list of permittivities}

  \item{polarisation}{polarisation}

  \item{...}{further args passed to multilayer}
}
\value{
  long format data.frame with positions and LFEF (para and
  perp)
}
\description{
  field profile in a ML stack
}
\details{
  runs multilayer and returns the LFIEF as a function of
  distance inside and outside of the structure
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
\seealso{
  Other helping_functions: \code{\link{invert_stack}}
}

