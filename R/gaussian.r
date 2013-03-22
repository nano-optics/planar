
##' Electric field from the transmission of a gaussian beam at a planar interface
##'
##' Integration is performed over a spectrum of incident plane waves using integrand_gb
##' @title gaussian_field
##' @param x position
##' @param y position
##' @param z position
##' @param lambda wavelength
##' @param alpha beam incident angle
##' @param psi beam polarisation angle
##' @param gamma beam waist radius
##' @param np incident medium (prism) refractive index
##' @param ns outer medium (substrate) refractive index
##' @return data.frame Electric field (squared modulus) at the xyz position
##' @export
##' @family gaussian_beam
##' @author Baptiste Auguie
gaussian_field <- function(x=1, y=1, z=1, 
                           lambda=500, alpha=pi/2 - 15*pi/180, psi=0, gamma=10, 
                           np=1.5, ns=1){

k0 <- 2*pi/lambda
kp <- k0*np

res <- adaptIntegrate(integrand_gb,
                      lowerLimit=c(-0.4, -0.4)*f, # some function of alpha
                      upperLimit=c(0.4, 0.4)*f, # some function of alpha
                      fDim = 6,
                      maxEval = 200,
                      r2 = c(x, y, z), kp=kp, psi=0, alpha=alpha,
                      gamma=gamma, np=np, ns=ns)$integral

E <- complex(real = res[1:3], imag=res[4:6])
Re(crossprod(E, Conj(E)))

}