
test_complex <- function(x){
  !isTRUE(all.equal(Im(x), 0))
}

##' epsilon_medium
##'
##' characterise the layers of a structure with unique labels for metals and dielectrics
##' @title epsilon_medium
##' @param epsilon list of real or complex values
##' @return factor
##' @export
##' @family user_level conversion utility
##' @author baptiste Auguie
epsilon_medium <- function(epsilon = list(3.5, 1, 3, 1, 12+1i, 3, 3.5)){
  
  is.metal <- sapply(epsilon, test_complex)
  metals <- unlist(unique(epsilon[is.metal]))
  dielectrics <- sort(unlist(unique(epsilon[!is.metal])))
  
  dnames <- paste0("n", seq_along(dielectrics))[match(epsilon, dielectrics)]
  mnames <- paste0("metal", seq_along(metals))[match(epsilon, metals)]
  
  factor(ifelse(is.na(dnames), mnames, dnames))
  
}

##' Local field intensity enhancement factors in a multilayer
##'
##' returns the LFIEFs as a function of distance inside and outside of the structure
##' @title lfief
##' @alias field_profile
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
lfief <- function(wavelength=500, angle=0, polarisation='p',
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
  m <- melt(all, id=1)
  
  ll = as.list(unique(m$L1))
  names(ll) = epsilon_medium(epsilon)
  m$material <- factor(m$L1)
  levels(m$material) = ll
  
  m
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
  positions <- c(-dmax, cumsum(thickness))
  positions[length(positions)] <- sum(thickness) + dmax
  n1 <- Re(sqrt(epsilon[1]))
  ## at least 10 points 
  ## stepsize <- min(c(diff(range(positions)) / res, 0.1*diff(positions)))
  stepsize <- diff(range(positions)) / res
  probes <- mapply(seq, positions[-length(positions)], positions[-1], 
                   MoreArgs=list(by=stepsize))
  id <- rep(seq_along(probes), sapply(probes, length))
  d <- unlist(probes)
  E <- planar$multilayer_field(k0, k0*sin(angle)*n1, unlist(epsilon),  
                          thickness, d, psi)$E
  if(field)
    return(E)
  
  ll = as.list(unique(id))
  names(ll) = epsilon_medium(epsilon)
  material <- factor(id)
  levels(material) = ll
  
  data.frame(x = d, I = Re(colSums(E*Conj(E))), id=id, material=material)
}
