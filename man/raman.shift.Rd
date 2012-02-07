\name{raman.shift}
\alias{raman.shift}
\title{raman.shift}
\usage{
  raman.shift(lambda = 500, shift = 1000, stokes = FALSE)
}
\arguments{
  \item{lambda}{wavelength (nm)}

  \item{shift}{Raman shift (cm-1)}

  \item{stokes}{logical Stokes or Anti-Stokes}
}
\value{
  wavelength of the Raman peak in nm
}
\description{
  Raman shift
}
\details{
  calculate the raman shift
}
\examples{
raman.shift(lambda=200)
}
\author{
  Baptiste Auguie
}

