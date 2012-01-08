\name{dipole.integrand}
\alias{dipole.integrand}
\title{dipole.integrand}
\usage{dipole.integrand(d=10, q, lambda, epsilon=list(incident = 1.5^2, 1^2),
    thickness=c(0, 0))}
\description{Dipole decay rates near a multilayer interface}
\details{Integrand of the dipole decay rates near a multilayer interface}
\author{baptiste Auguie}
\arguments{\item{d}{distance in nm}
\item{q}{normalised in-plane wavevector in [0, infty) OR sqrt(1-q^2) if change.variable}
\item{lambda}{wavelength in nm}
\item{epsilon}{list of dielectric functions}
\item{thickness}{list of layer thicknesses}}
