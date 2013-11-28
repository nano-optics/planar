
##' invert the description of a multilayer to simulate the opposite direction of incidence
##'
##' inverts list of epsilon and thickness of layers
##' @title invert_stack
##' @param p list
##' @return list
##' @export
##' @family helping_functions
##' @author Baptiste Auguie
invert_stack <- function(p){
  p[["epsilon"]] <- rev(p[["epsilon"]])
  p[["thickness"]] <- rev(p[["thickness"]])
  p
}

##' relabel factors
##'
##' Wide to long format data.frame with new factor variable(s) describing the original columns
##' @title classify
##' @param d data.frame
##' @param id column id
##' @param vars variables
##' @param ... passed on to melt
##' @return data.frame
##' @export
##' @family helping_functions
##' @author Baptiste Auguie
classify <- function(d, id=NULL, vars=NULL, ...){

  m <- melt(d, id.vars=id, ...)

  id.variables <- list()
  for (ii in seq_along(vars)){
    id.variables[[ii]] <- rep(vars[[ii]], each=nrow(d))
  }
  names(id.variables) <- names(vars)

  data.frame(m, id.variables)
}

##' relabel factors
##'
##' @title modify_levels
##' @param f factor
##' @param modify named list
##' @return factor
##' @export
##' @family helping_functions
##' @author Baptiste Auguie
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
##' @param stokes logical Stokes or Anti-Stokes
##' @return wavelength of the Raman peak in nm
##' @author Baptiste Auguie
raman_shift <- function(wavelength = 632.8, shift = 520, stokes = TRUE){

  if(stokes) 1 / (1/wavelength - shift * 1e-7) else 1 / (1/wavelength + shift * 1e-7)

}
