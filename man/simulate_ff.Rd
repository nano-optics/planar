\name{simulate_ff}
\alias{simulate_ff}
\title{simulate_ff}
\usage{
simulate_ff(..., s = NULL, fun = tamm_stack, wavelength = seq(400, 1000),
  angle = 0, polarisation = c("p", "s"))
}
\arguments{
  \item{fun}{function returning a stack}

  \item{wavelength}{numeric vector}

  \item{angle}{incident angle in radians}

  \item{polarisation}{p or s}
}
\value{
data.frame
}
\description{
simultate the far-field response of a multilayer stack
}
\details{
wrapper around recursive_fresnelcpp for a stack structure
}
\author{
Baptiste Auguie
}
\seealso{
Other helping_functions user_level stack:
\code{\link{rev.stack}}; \code{\link{simulate_nf}}
}

