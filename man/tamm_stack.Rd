\name{tamm_stack}
\alias{tamm_stack}
\title{tamm_stack}
\usage{
tamm_stack(lambda0 = 630, n1 = 1.28, n2 = 1.72, d1 = lambda0/4/n1,
  d2 = lambda0/4/n2, N = 2 * pairs, pairs = 4, dx1 = 0, dx2 = 0,
  nx1 = n1, nx2 = n2, dm = 50, metal = "epsAu", position = c("after",
  "before"), nleft = 1.5, nright = 1, ...)
}
\arguments{
  \item{lambda0}{central wavelength of the stopband}

  \item{n1}{real refractive index for odd layers}

  \item{n2}{real refractive index for even layers}

  \item{d1}{odd layer thickness in nm}

  \item{d2}{even layer thickness in nm}

  \item{N}{number of layers, overwrites pairs}

  \item{pairs}{number of pairs}

  \item{nx1}{real refractive index for last odd layer}

  \item{nx2}{real refractive index for last even layer}

  \item{dx1}{variation of last odd layer thickness in nm}

  \item{dx2}{variation of last even layer thickness in nm}

  \item{dm}{thickness of metal layer}

  \item{metal}{character name of dielectric function}

  \item{nleft}{refractive index of entering medium}

  \item{nright}{refractive index of outer medium}

  \item{position}{metal position relative to DBR}

  \item{...}{ignored}
}
\value{
list of class 'stack'
}
\description{
DBR-metal stack structure
}
\details{
periodic structure of dielectric layers against metal film
}
\author{
baptiste Auguie
}
\seealso{
Other stack user_level: \code{\link{dbr_stack}},
\code{\link{embed_stack}}, \code{\link{layer_stack}}
}

