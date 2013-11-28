
##' Electric field from the transmission of a gaussian beam at a planar interface
##'
##' Integration is performed over a spectrum of incident plane waves using integrand_gb
##' @title gaussian_field
##' @param x position
##' @param y position
##' @param z position
##' @param wavelength wavelength
##' @param alpha beam incident angle
##' @param psi beam polarisation angle
##' @param w0 beam waist radius
##' @param ni incident medium (prism) refractive index
##' @param no outer medium (substrate) refractive index
##' @param cutoff radial integration limit
##' @param maxEval passed to adaptIntegrate
##' @return data.frame Electric field (squared modulus) at the x, y, z position
##' @export
##' @family gaussian_beam
##' @author Baptiste Auguie
gaussian_field <- function(x=1, y=1, z=1, wavelength=500, alpha = 15*pi/180, psi=0, 
                           w0=1e4, ni=1.5, no=1,
                           cutoff = min(1, sqrt(3*4)/(w0*Re(2*pi/wavelength*ni))),
                           maxEval = 300){
  
  k <- 2*pi/wavelength
  ki <- k*ni
  res <- adaptIntegrate(gaussian$integrand_gb,
                        lowerLimit=c(0, 0), # rho in [0,1], angle in [0,2*pi]
                        upperLimit=c(cutoff, 2*pi), # exp(-3) is tiny
                        fDim = 6, tol = 1e-04,
                        maxEval = maxEval,
                        r2 = c(x, y, z), ki=ki, psi=psi, alpha=alpha,
                        w0=w0, ni=ni, no=no)$integral
  
  E <- complex(real = res[1:3], imag=res[4:6])
  Re(crossprod(E, Conj(E)))
  
}



