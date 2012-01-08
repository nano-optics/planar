\name{field.profile}
\alias{field.profile}
\title{field.profile}
\usage{field.profile(lambda=500, theta=0, polarisation="p", thickness=c(0, 20,
    140, 20, 0), dmax=200, res=1000, epsilon=list(1^2, -12 + (0 + (0+1i)), 1.38^2, -12 + (0 + (0+1i)), 1.46^2), ...)}
\description{field profile in a ML stack}
\details{runs multilayer and returns the LFEF as a function of distance inside and outside of the structure}
\value{long format data.frame with positions and LFEF (para and perp)}
\author{baptiste Auguie}
\arguments{\item{lambda}{wavelength}
\item{theta}{angle}
\item{dmax}{maximum distance to interface, if > layer thickness}
\item{thickness}{vector of layer thickness}
\item{res}{resolution of sampling points}
\item{epsilon}{list of permittivities}
\item{polarisation}{polarisation}
\item{...}{further args passed to multilayer}}
