
integration_points_disk <- function(Nquad = 30, 
                                    min=c(0, 0), max=c(1, 2*pi),
                                    method="GL",  init=TRUE){
  
  method <- match.arg(method)
  rho1=min[1]; rho2=max[1]; theta1=min[2]; theta2=max[2]; 
  C1 = (rho2 - rho1) / 2;  D1 = (rho1 + rho2) / 2;
  C2 = (theta2 - theta1) / 2; D2 = (theta1 + theta2) / 2;
  
  rndN <- ceiling(sqrt(Nquad))
  GL_rho <- statmod::gauss.quad(rndN)
  GL_theta <- statmod::gauss.quad(rndN)
  
  rho = GL_rho$nodes*C1 + D1  
  theta = GL_theta$nodes*C2 + D2
  
  # grid of params
  grid <- expand.grid(rho=rho, theta=theta)
  # corresponding weights for 2D quadrature
  weights <- expand.grid(rho=GL_rho$weights, theta=GL_theta$weights)
  # combine the weigths for each point; grid$rho comes from the Jacobian in the integral
  weights <- C1 * C2 * grid$rho * weights$rho * weights$theta
  
  return(list(grid=as.matrix(grid), weights=weights))
}  


##' Electric field from the transmission of a gaussian beam at a planar interface
##'
##' Integration is performed over a spectrum of incident plane waves using integrand_gb
##' @title gaussian_near_field
##' @param x position
##' @param y position
##' @param z position
##' @param wavelength wavelength
##' @param alpha beam incident angle
##' @param psi beam polarisation angle
##' @param w0 beam waist radius
##' @param ni incident medium (prism) refractive index
##' @param nl layer refractive index
##' @param no outer medium (substrate) refractive index
##' @param d thickness of layer
##' @param cutoff radial integration limit
##' @param maxEval passed to adaptIntegrate
##' @param tol passed to adaptIntegrate
##' @param field logical: return the electric field (complex vector), or modulus squared
##' @return data.frame electric field at the x, y, z position
##' @export
##' @family gaussian_beam
##' @author Baptiste Auguie
gaussian_near_field <- function(x=1, y=1, z=1, wavelength=500, alpha = 15*pi/180, psi=0, 
                           w0=1e4, ni=1.5, no=1.0, nl=no, d=0,
                           cutoff = min(1, 3*wavelength/(ni*pi*w0)), # integrand ~ exp(-3^2)
                           maxEval = 3000, tol=1e-04, field=FALSE){
  
  k <- 2*pi/wavelength
  ki <- k*ni
  res <- cubature::adaptIntegrate(gaussian$integrand_gb,
                                  lowerLimit=c(0, 0), # rho in [0,1], angle in [0,2*pi]
                                  upperLimit=c(cutoff, 2*pi), 
                                  fDim = 6, tol = tol,
                                  maxEval = maxEval,
                                  r2 = c(x, y, z), ki=ki, psi= psi, alpha=alpha,
                                  w0=w0, ni=ni, no=no, nl=nl, d=d)$integral
  
  E <- complex(real = res[1:3], imaginary=res[4:6])
  if(field) return(E)
  
  Re(crossprod(E, Conj(E)))
}




##' Electric field from the transmission of a gaussian beam at a planar interface
##'
##' Integration is performed over a spectrum of incident plane waves using integrand_gb
##' @title gaussian_near_field2
##' @param x position
##' @param y position
##' @param z position
##' @param wavelength wavelength
##' @param alpha beam incident angle
##' @param psi beam polarisation angle
##' @param w0 beam waist radius
##' @param epsilon vector of permittivities
##' @param thickness thickness corresponding to each medium
##' @param cutoff radial integration limit
##' @param maxEval passed to adaptIntegrate
##' @param tol passed to adaptIntegrate
##' @param field logical: return the electric field (complex vector), or modulus squared
##' @return data.frame electric field at the x, y, z position
##' @export
##' @family gaussian_beam
##' @author Baptiste Auguie
gaussian_near_field2 <- function(x=1, y=1, z=1, wavelength=632.8, alpha = 15*pi/180, psi=0, 
                                 w0=1e4, epsilon = c(1.5^2, epsAg(lambda)$epsilon, 1.0^2, 1.0^2),
                                 thickness = c(0, 50, 10, 0),
                                 cutoff = min(1, 3*wavelength/(Re(sqrt(epsilon[1]))*pi*w0)), 
                                 maxEval = 3000, tol=1e-04, field=FALSE){
  
  k0 <- 2*pi/wavelength
  res <- cubature::adaptIntegrate(gaussian$integrand_gb2,
                                  lowerLimit=c(0, 0), # rho in [0,1], angle in [0,2*pi]
                                  upperLimit=c(cutoff, 2*pi), 
                                  fDim = 6, tol = tol,
                                  maxEval = maxEval,
                                  r2 = c(x, y, z), k0=k0, psi= psi, alpha=alpha,
                                  w0=w0, epsilon=epsilon, thickness=thickness)$integral
  
  E <- complex(real = res[1:3], imaginary=res[4:6])
  if(field) return(E)
  
  Re(crossprod(E, Conj(E)))
}

