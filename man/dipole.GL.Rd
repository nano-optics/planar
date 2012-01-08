\name{dipole.GL}
\alias{dipole.GL}
\title{dipole2}
\usage{dipole.GL(d=1, lambda, epsilon=list(incident = 1^2), thickness=c(0, 0),
    qcut, show.messages=FALSE, Nquadrature1=50, Nquadrature2=200,
    Nquadrature3=50)}
\description{Dipole decay rates near a multilayer interface}
\details{dipole decay rates near a multilayer interface}
\author{baptiste Auguie}
\arguments{\item{d}{distance in nm}
\item{lambda}{wavelength in nm}
\item{epsilon}{list of dielectric functions}
\item{thickness}{list of layer thicknesses}
\item{Nquadrature1}{quadrature points in radiative region}
\item{Nquadrature2}{quadrature points in SPPs region}
\item{Nquadrature3}{quadrature points in dipole image region}
\item{qcut}{transition between regions 2 and 3}
\item{show.messages}{logical, display progress bar}}
