\name{integrand3}
\alias{integrand3}
\title{integrand3}
\usage{integrand3(d=10, u, ucut, lambda, epsilon=list(incident = 1.5^2, 1^2),
    thickness=c(0, 0))}
\description{Dipole decay rates near a multilayer interface}
\details{Integrand of the dipole decay rates near a multilayer interface. Transformed part III
from u=ucut to infinity}
\author{baptiste Auguie}
\arguments{\item{d}{distance in nm}
\item{u}{transformed normalised in-plane wavevector sqrt(q^2 - 1)}
\item{ucut}{limit of the integral}
\item{lambda}{wavelength in nm}
\item{epsilon}{list of dielectric functions}
\item{thickness}{list of layer thicknesses}}
