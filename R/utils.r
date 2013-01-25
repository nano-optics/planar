
## Wide to long format data.frame with new factor variable(s) describing the original columns
classify <- function(d, id=NULL, vars=NULL, ...){

  m <- melt(d, id.vars=id, ...)

  id.variables <- list()
  for (ii in seq_along(vars)){
    id.variables[[ii]] <- rep(vars[[ii]], each=nrow(d))
  }
  names(id.variables) <- names(vars)

  data.frame(m, id.variables)
}

modify_levels <- function(f, modify=list()){
  f <- factor(f)
  levs = levels(f)
  m = match(modify,levs)
  levs[m] = names(modify)
  factor(f,labels=levs)
}


Curry <- function (FUN, ...) 
{
    .orig = list(...)
    function(...) do.call(FUN, c(.orig, list(...)))
}

##' Raman shift
##'
##' Raman shift conversion to absolute wavelength
##' @title raman_shift
##' @export
##' @param wavelength wavelength (nm)
##' @param shift Raman shift (cm-1)
##' @param antistokes logical Stokes or Anti-Stokes
##' @return wavelength of the Raman peak in nm
##' @author Baptiste Auguie
##' @examples
##' raman.shift(wavelength=200)
raman_shift <- function(wavelength = 632.8, shift = 520, stokes = TRUE){

  if(stokes) 1 / (1/wavelength - shift * 1e-7) else 1 / (1/wavelength + shift * 1e-7)

}
