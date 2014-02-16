\name{simulate_nf}
\alias{simulate_nf}
\title{simulate_nf}
\usage{
simulate_nf(..., s = NULL, fun = tamm_stack, wavelength = 630,
  angle = 0, polarisation = c("p", "s"), res = 10000, field = FALSE)
}
\arguments{
  \item{fun}{function returning a stack}

  \item{wavelength}{numeric vector}

  \item{angle}{incident angle in radians}

  \item{polarisation}{p or s}

  \item{res}{number of points}

  \item{field}{logical, return the complex electric field
  vector}
}
\value{
data.frame
}
\description{
simultate the internal field of a multilayer stack
}
\details{
wrapper around multilayer_field for a stack structure
}
\author{
Baptiste Auguie
}
\seealso{
Other helping_functions user_level stack:
\code{\link{rev.stack}}; \code{\link{simulate_ff}}
}

