% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stacks.r
\name{dbr_stack}
\alias{dbr_stack}
\title{dbr_stack}
\usage{
dbr_stack(
  lambda0 = 630,
  n1 = 1.28,
  n2 = 1.72,
  d1 = lambda0/4/n1,
  d2 = lambda0/4/n2,
  N = 2 * pairs,
  pairs = 4,
  ...
)
}
\arguments{
\item{lambda0}{central wavelength of the stopband}

\item{n1}{real refractive index for odd layers}

\item{n2}{real refractive index for even layers}

\item{d1}{odd layer thickness in nm}

\item{d2}{even layer thickness in nm}

\item{N}{number of layers, overwrites pairs}

\item{pairs}{number of pairs}

\item{...}{ignored}
}
\value{
list of class 'stack'
}
\description{
DBR stack structure
}
\details{
periodic structure of dielectric layers
}
\seealso{
Other stack user_level: 
\code{\link{embed_stack}()},
\code{\link{layer_stack}()},
\code{\link{tamm_stack_ir}()},
\code{\link{tamm_stack_porous}()},
\code{\link{tamm_stack}()}
}
\author{
baptiste Auguie
}
\concept{stack user_level}
