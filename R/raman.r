
##' Raman shift
##'
##' calculate the raman shift
##' @title raman.shift
##' @export
##' @param lambda wavelength (nm)
##' @param shift Raman shift (cm-1)
##' @param stokes logical Stokes or Anti-Stokes
##' @return wavelength of the Raman peak in nm
##' @author Baptiste Auguie
##' @examples
##' raman.shift(lambda=200)
raman.shift <- function(lambda=500, shift=1000, stokes=FALSE){

  lambdaR <- 1 / shift * 1e7

  if(stokes) 1 / (1/lambda + 1/lambdaR) else 1 / (1/lambda - 1/lambdaR)

}
