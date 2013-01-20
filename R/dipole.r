## Calculates Mrad and Mtot for a dipole near a multilayer, using the angular decomposition of the dipole field into plane waves
## Mrad is done in one step, but the integration for Mtot is divided in 3 regions
## thus we define 4 different integrands


##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the radiative dipole decay rates near a multilayer interface. 
##' @title integrand_rad
##' @export
##' @param d distance in nm
##' @param theta angle in radians
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @author baptiste Auguie
integrand_rad <- function(d = 10, theta, lambda,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0)){

  ## for 0 < q < 1, i.e. 0 < u < 1
  
  ## define constants
  k0 <- 2*pi/lambda
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Ntheta <- length(theta)
  
  cost <- cos(theta)
  sint <- sin(theta)
  
  rp <- recursive_fresnelcpp(lambda=lambda,
                           q = sint,
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(lambda=lambda,
                           q = sint,
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection

  phase <- exp(2i*d*outer(k1,cost))
  
  integrand.p <- Mod(matrix(1, Nlambda, Ntheta, byrow=TRUE) + rp * phase)^2 *
    matrix(sint^3, Nlambda, Ntheta, byrow=TRUE)
  
  integrand.s <- (Mod(matrix(1, Nlambda, Ntheta, byrow=TRUE) + rs * phase)^2 +
                  Mod(matrix(1, Nlambda, Ntheta, byrow=TRUE) - rp * phase)^2 *
                  matrix(cost^2, Nlambda, Ntheta, byrow=TRUE)) *
                    matrix(sint, Nlambda, Ntheta, byrow=TRUE)
      
  rbind(integrand.p, integrand.s)
}


##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface. Transformed part I (radiative)
##' from u=0 to 1
##' @title integrand_nr1
##' @export
##' @param d distance in nm
##' @param u transformed normalised in-plane wavevector sqrt(1-q^2)
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @author baptiste Auguie
integrand_nr1 <- function(d=10, u, lambda,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0)){

  ## integrand1 is for 0 < q < 1, i.e. 0 < u < 1
  
  ## define constants
  k0 <- 2*pi/lambda
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(u)
  
  rp <- recursive_fresnelcpp(lambda=lambda,
                           q = sqrt(1 - u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(lambda=lambda,
                           q = sqrt(1 - u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection

    phase <- exp(2i*d*outer(k1,u))
    ## print(c(Nlambda, Nq, dim(rp), dim(phase)))
    integrand.p <- Re(matrix(1 - u^2, Nlambda, Nq, byrow=TRUE) * rp * phase)
    integrand.s <- Re(( rs - rp*matrix(u^2, Nlambda, Nq, byrow=TRUE)) * phase)
      
  rbind(integrand.p, integrand.s)
}

##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface. Transformed part II
##' from u=0 to ucut
##' @title integrand_nr2
##' @export
##' @param d distance in nm
##' @param u transformed normalised in-plane wavevector sqrt(q^2 - 1)
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @author baptiste Auguie
integrand_nr2 <- function(d=10, u, lambda,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0)){

  ## integrand2 is for 1 < q < infty, i.e. 0 < u < infty

  ## define constants
  k0 <- 2*pi/lambda
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(u)
  
  rp <- recursive_fresnelcpp(lambda=lambda,
                           q = sqrt(1 + u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(lambda=lambda,
                           q = sqrt(1 + u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection
  
  ## phase is now purely real
  phase <- exp(-2*d*outer(k1,u))

  ## take the imaginary part now
  integrand.p <- matrix(1 + u^2, Nlambda, Nq, byrow=TRUE) * Im( rp ) * phase
  integrand.s <- Im(rs + rp*matrix(u^2, Nlambda, Nq, byrow=TRUE)) * phase
  
  rbind(integrand.p, integrand.s)
}


##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface. Transformed part III
##' from u=ucut to infinity
##' @title integrand_nr3
##' @export
##' @param d distance in nm
##' @param u transformed normalised in-plane wavevector sqrt(q^2 - 1)
##' @param ucut limit of the integral
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @author baptiste Auguie
integrand_nr3 <- function(d=10, u, ucut, lambda,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0)){

  ## define constants
  k0 <- 2*pi/lambda
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(u)
  
  ## integrand2 is for ucut < u < infty
  ## performing a change of variable mapping u in [ucut, infty) -> [0,1]
  ## change of variables
  ## \int_a^\infty f(x)dx = \int_0^1 f(a + t/(1-t)). 1 / (1-t)^2 dt
  ## as suggested on http://ab-initio.mit.edu/wiki/index.php/Cubature

  ## new variable
  t <-  ucut + u / (1 - u)
  ## Jacobian of transformation
  Jac <-  matrix(1 / (1 - u)^2, Nlambda, Nq, byrow=TRUE)
  
  
  rp <- recursive_fresnelcpp(lambda=lambda,
                           q = sqrt(1 + t^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(lambda=lambda,
                           q = sqrt(1 + t^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection
  
  ## phase is now purely real
  phase <- exp(-2*d*outer(k1,t))

  ## take the imaginary part now
  integrand.p <- matrix(1 + t^2, Nlambda, Nq, byrow=TRUE) * Im( rp ) * phase * Jac
  integrand.s <- Im(rs + rp*matrix(t^2, Nlambda, Nq, byrow=TRUE)) * phase * Jac
  
  rbind(integrand.p, integrand.s)
}

##' Dipole decay rates near a multilayer interface
##'
##' dipole decay rates near a multilayer interface
##' @title dipole2
##' @export
##' @param d distance in nm
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @param Nquadrature1 maximum number of quadrature points in radiative region
##' @param Nquadrature2 maximum number of quadrature points in SPPs region
##' @param Nquadrature3 maximum number of quadrature points in dipole image region
##' @param rel.err relative error
##' @param show.messages logical, display progress bar
##' @param qcut transition between regions 2 and 3
##' @author baptiste Auguie
dipole <- function(d=1,
                   lambda,
                   epsilon = list(incident=1.0^2),
                   thickness = c(0, 0), qcut=NULL, rel.err = 1e-3,
                   Nquadrature1 = 1e3, Nquadrature2 = 1e4, Nquadrature3 = 1e4,
                   show.messages=TRUE){
   
  Nlambda <- length(lambda)
  
  require(cubature) # quadrature

  ## if no qcut provided, estimate one from max of
  ## all possible SPP dispersions
  if(is.null(qcut)){
    qcut <- 1.1

    epsilon_norm <- do.call(cbind, epsilon)
    
    for(ii in seq(1, length(epsilon) - 1)){
      qspp <- sqrt(epsilon_norm[,ii] / epsilon_norm[,1])*
        sqrt(epsilon_norm[,ii+1] / (epsilon_norm[,ii] + epsilon_norm[,ii+1]))
      
      qcut <- max(qcut, max(Re(qspp)))
    }

    if(show.messages)
      print(paste("using qcut=",round(qcut,2)))
    
  }
  
  ## integration from 0 to 1 for the transformed radiative bit
  
  in1 <- adaptIntegrate(integrand_nr1, lowerLimit = 0,
                        upperLimit = 1,
                        fDim = 2*Nlambda, tol=rel.err,
                        maxEval = Nquadrature1, 
                        d=d, lambda=lambda,
                        epsilon=epsilon, thickness=thickness)

  integral1.perp <- in1$integral[seq(1,Nlambda)]
  integral1.par <- in1$integral[seq(Nlambda+1,2*Nlambda)]
  
  ## integration from 0 to ucut
  
  ucut <- sqrt(qcut^2 - 1)
  
  in2 <- adaptIntegrate(integrand_nr2, lowerLimit = 0,
                        upperLimit = ucut,
                        fDim = 2*Nlambda, tol=rel.err,
                        maxEval = Nquadrature2,
                        d=d, lambda=lambda,
                        epsilon=epsilon, thickness=thickness)

  integral2.perp <- in2$integral[seq(1,Nlambda)]
  integral2.par <- in2$integral[seq(Nlambda+1,2*Nlambda)]
  
  ## integration from ucut to Inf
  ## integrand performing a change of variable mapping u in [ucut, infty) -> t in [0,1]
  
  in3 <- adaptIntegrate(integrand_nr3, lowerLimit = 0,
                        upperLimit = 1,
                        fDim = 2*Nlambda, tol=rel.err,
                        maxEval = Nquadrature3, 
                        ucut=ucut, d=d, lambda=lambda,
                        epsilon=epsilon, thickness=thickness)

  integral3.perp <- in3$integral[seq(1,Nlambda)]
  integral3.par <- in3$integral[seq(Nlambda+1,2*Nlambda)]
  
  ## Mrad
  
  in4 <- adaptIntegrate(integrand_rad, lowerLimit = 0,
                        upperLimit = pi/2,
                        fDim = 2*Nlambda, tol=rel.err,
                        maxEval = Nquadrature1, 
                        d=d, lambda=lambda,
                        epsilon=epsilon, thickness=thickness)
  
  Mrad.perp <- in4$integral[seq(1,Nlambda)]
  Mrad.par <- in4$integral[seq(Nlambda+1,2*Nlambda)]

  evaluations <- c(in1$functionEvaluations,
                   in2$functionEvaluations,
                   in3$functionEvaluations,
                   in4$functionEvaluations)
  if(show.messages)
    message(c("integration points=", paste(round(evaluations,2))))
  
  list(results= 
       data.frame(wavelength=lambda,
                  Mtot.perp = 1 + 3/2*(integral1.perp + integral2.perp + integral3.perp),
                  Mtot.par = 1 + 3/4*(integral1.par + integral2.par + integral3.par),
                  Mrad.perp = 3/4 * Mrad.perp, Mrad.par = 3/8 * Mrad.par) ,
       errors = list(in1$error, in2$error, in3$error, in4$error),
       evaluations = evaluations)

}

