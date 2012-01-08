\name{integrand2}
\alias{integrand2}
\title{integrand2}
\usage{integrand2(d=10, u, lambda, epsilon=list(incident = 1.5^2, 1^2),
    thickness=c(0, 0))}
\description{Dipole decay rates near a multilayer interface}
\details{Integrand of the dipole decay rates near a multilayer interface. Transformed part II
from u=0 to ucut}
\author{baptiste Auguie}
\arguments{\item{d}{distance in nm}
\item{u}{transformed normalised in-plane wavevector sqrt(q^2 - 1)}
\item{lambda}{wavelength in nm}
\item{epsilon}{list of dielectric functions}
\item{thickness}{list of layer thicknesses}}
