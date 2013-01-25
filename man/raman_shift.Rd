\name{raman_shift}
\alias{raman_shift}
\title{raman_shift}
\usage{
  raman_shift(wavelength = 632.8, shift = 520,
    stokes = TRUE)
}
\arguments{
  \item{wavelength}{wavelength (nm)}

  \item{shift}{Raman shift (cm-1)}

  \item{antistokes}{logical Stokes or Anti-Stokes}
}
\value{
  wavelength of the Raman peak in nm
}
\description{
  Raman shift
}
\details{
  Raman shift conversion to absolute wavelength
}
\examples{
raman.shift(wavelength=200)
}
\author{
  Baptiste Auguie
}

