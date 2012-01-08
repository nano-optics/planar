\name{multilayer2}
\alias{multilayer2}
\title{multilayer2}
\usage{multilayer2(lambda=NULL, k0=2 * pi/lambda, theta=NULL, q=sin(theta),
    epsilon=list(incident = 1.5^2, 1.33), thickness=c(0, 0),
    polarisation=c("p", "s"))}
\description{Multilayer Fresnel coefficients}
\details{solves the EM problem of a multilayered interface}
\value{fresnel coefficients and field profiles}
\author{baptiste Auguie}
\references{Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).}
\arguments{\item{lambda}{[vector] wavelength in nm}
\item{k0}{[vector] wavevector in nm^-1}
\item{theta}{[vector] incident angles in radians}
\item{q}{[vector] normalised incident in-plane wavevector}
\item{epsilon}{list of N+2 dielectric functions, each of length 1 or length(lambda)}
\item{thickness}{vector of N+2 layer thicknesses, first and last are dummy}
\item{polarisation}{[character] switch between p- and s- polarisation}}
\examples{library(planar)
demo(package="planar")}
