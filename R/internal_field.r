##' field profile in a ML stack
##'
##' returns the electric field as a function of distance inside and outside of the structure
##' @title field_profile
##' @export
##' @param wavelength wavelength
##' @param angle angle
##' @param dmax maximum distance to interface, if > layer thickness
##' @param thickness vector of layer thickness
##' @param res resolution of sampling points
##' @param res2 resolution of sampling points outside stack
##' @param epsilon list of permittivities
##' @param polarisation polarisation
##' @param displacement logical, Mperp corresponds to displacement squared (D=epsilon x E)
##' @param ... further args passed to multilayer 
##' @return long format data.frame with positions and LFEF (para and perp)
##' @author baptiste Auguie
##' @family helping_functions
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
##' 
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
field_profile <- function(wavelength=500, angle=0, polarisation='p',
                          thickness = c(0, 20, 140, 20, 0), 
                          dmax=200,  res=1e3, res2=res/10,
                          epsilon=list(1^2, -12 , 1.38^2, -12 , 1.46^2), 
                          displacement=FALSE, ...){
  
  d <- seq(0, max(thickness), length=res)
  dout <- seq(0, dmax, length=res2)
  res <- multilayer(wavelength=wavelength, angle=angle,
                    epsilon=epsilon,
                    thickness = thickness, d=d, dout=dout,
                    polarisation=polarisation, ...)
  
  #   loop over 2nd to last interface
  all <- lapply(seq(1,length(res$dist) - 1), function(lay){
    Mperp <- if(displacement) (Mod(epsilon[[lay+1]]))^2 * res$Mr.perp[[lay]] else res$Mr.perp[[lay]]
    data.frame(x = res$dist[[lay+1]], M.par=res$Mr.par[[lay]],
               M.perp=Mperp)
  })
  #   combine with first interface
  Mperp <- if(displacement) (Mod(epsilon[[1]]))^2 * res$Ml.perp[[1]] else res$Ml.perp[[1]]
  all <- c(list(data.frame(x = res$dist[[1]],
                           M.par=res$Ml.par[[1]],
                           M.perp=Mperp)), all)
  
  names(all) <- paste("layer", seq_along(res$dist))
  melt(all, id=1)
}

##' Internal field in a ML stack
##'
##' returns the electric field as a function of distance inside and outside of the structure
##' @title internal_field
##' @export
##' @param wavelength wavelength
##' @param angle angle
##' @param psi polarisation angle (0 for TM)
##' @param dmax maximum distance to interface
##' @param thickness vector of layer thickness
##' @param res resolution of sampling points
##' @param epsilon permittivities
##' @param field logical, return complex field vector, or modulus squared
##' @param ... further args ignored
##' @return data.frame with position and electric field vector
##' @author baptiste Auguie
##' @family helping_functions
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
##' 
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
internal_field <- function(wavelength=500, angle=0, psi=0,
                           thickness = c(0, 20, 140, 20, 0), 
                           dmax=200,  res=1e3, 
                           epsilon=c(1^2, -12 , 1.38^2, -12 , 1.46^2), 
                           field = FALSE, ...){
  k0 <- 2*pi/wavelength
  positions <- c(-dmax, thickness)
  positions[length(positions)] <- cumsum(thickness) + dmax
  n1 <- Re(sqrt(epsilon[1]))
  d <- unique(c(mapply(seq, positions[-length(positions)], positions[-1], length=res)))
  E <- planar$multilayer_field(k0, k0*sin(angle)*n1, unlist(epsilon),  
                          thickness, d, psi)
  if(field)
    return(E)
  
  data.frame(x = d, I = Re(colSums(E*Conj(E))))
}
