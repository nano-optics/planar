
##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the radiative dipole decay rates near a multilayer interface
##' @title integrand.radGL
##' @export
##' @param d distance in nm
##' @param theta angle in radians
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @author baptiste Auguie
integrand.radGL <- function(d=10, theta, lambda,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0)){

  ## integrand1 is for 0 < q < 1, i.e. 0 < u < 1
  
  ## define constants
  k0 <- 2*pi/lambda
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Ntheta <- length(theta)
  
  cost <- cos(theta)
  sint <- sin(theta)
  
  rp <- recursive.fresnel2(lambda=lambda,
                           q = sint,
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive.fresnel2(lambda=lambda,
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
      
  list(integrand.p = integrand.p, integrand.s = integrand.s)
}


##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface
##' @title integrand1GL
##' @export
##' @param d distance in nm
##' @param u transformed normalised in-plane wavevector sqrt(1-q^2)
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @author baptiste Auguie
integrand1GL <- function(d=10, u, lambda,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0)){

  ## integrand1 is for 0 < q < 1, i.e. 0 < u < 1
  
  ## define constants
  k0 <- 2*pi/lambda
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(u)
  
  rp <- recursive.fresnel2(lambda=lambda,
                           q = sqrt(1 - u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive.fresnel2(lambda=lambda,
                           q = sqrt(1 - u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection

    phase <- exp(2i*d*outer(k1,u))
    ## print(c(Nlambda, Nq, dim(rp), dim(phase)))
    integrand.p <- Re(matrix(1 - u^2, Nlambda, Nq, byrow=TRUE) * rp * phase)
    integrand.s <- Re(( rs - rp*matrix(u^2, Nlambda, Nq, byrow=TRUE)) * phase)
      
  list(integrand.p = integrand.p, integrand.s = integrand.s)
}


##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface
##' @title integrand2GL
##' @export
##' @param d distance in nm
##' @param u transformed normalised in-plane wavevector sqrt(q^2 - 1)
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @author baptiste Auguie
integrand2GL <- function(d=10, u, lambda,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0)){

  ## integrand2 is for 1 < q < infty, i.e. 0 < u < infty
  
  ## define constants
  k0 <- 2*pi/lambda
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(u)
  
  rp <- recursive.fresnel2(lambda=lambda,
                           q = sqrt(1 + u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive.fresnel2(lambda=lambda,
                           q = sqrt(1 + u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection
  
  ## phase is now purely real
  phase <- exp(-2*d*outer(k1,u))

  ## take the imaginary part now
  integrand.p <- matrix(1 + u^2, Nlambda, Nq, byrow=TRUE) * Im( rp ) * phase
  integrand.s <- Im(rs + rp*matrix(u^2, Nlambda, Nq, byrow=TRUE)) * phase
  
  list(integrand.p = integrand.p, integrand.s = integrand.s)
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
##' @param Nquadrature1 quadrature points in radiative region
##' @param Nquadrature2 quadrature points in SPPs region
##' @param Nquadrature3 quadrature points in dipole image region
##' @param qcut transition between regions 2 and 3
##' @param show.messages logical, display progress bar
##' @author baptiste Auguie
dipole.GL <- function(d=1,
                   lambda,
                   epsilon = list(incident=1.0^2),
                   thickness = c(0, 0), qcut=NULL, show.messages=FALSE,
                   Nquadrature1 = 50, Nquadrature2 = 200, Nquadrature3 = 50){
   
  Nlambda <- length(lambda)
  
  require(statmod) # quadrature points in (-1, 1)

  GL1 <- gauss.quad(Nquadrature1)
  GL2 <- gauss.quad(Nquadrature2)
  GL3 <- gauss.quad(Nquadrature3)

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
  umax1 <- 1; umin1 <- 0;
  C1 <- (umax1 - umin1)/2 ; D1 <- (umax1+umin1)/2
  unodes1 <- C1 * GL1$nodes + D1
  uweights1 <- GL1$weights * C1

  Nu1 <- length(unodes1)
  
  in1 <- integrand1GL(u=unodes1,
                    d=d, lambda=lambda,
                    epsilon=epsilon, thickness=thickness)
    
  weights1 <- matrix(uweights1, nrow=Nlambda, ncol=Nu1, byrow=TRUE)

  integral1.perp <- rowSums(in1$integrand.p*weights1)
  integral1.par <- rowSums(in1$integrand.s*weights1)
  
  ## integration from 1 to ucut
  
  ucut <- sqrt(qcut^2 - 1)
  
  umax2 <- ucut; umin2 <- 0;
  C2 <- (umax2 - umin2)/2 ; D2 <- (umax2 + umin2)/2
    
  unodes2 <- C2 * GL2$nodes + D2
  uweights2 <- GL2$weights * C2 
    
  Nu2 <- length(unodes2)
  
  in2 <- integrand2GL(u=unodes2,
                    d=d, lambda=lambda,
                    epsilon=epsilon, thickness=thickness)
  
  weights2 <- matrix(uweights2, nrow=Nlambda, ncol=Nu2, byrow=TRUE)
  
  integral2.perp <- rowSums(in2$integrand.p*weights2)
  integral2.par <- rowSums(in2$integrand.s*weights2)
  
  ## integration from ucut to Inf
  ## performing a change of variable mapping u in [ucut, infty) -> [0,1]
  ## change of variables
  ## \int_a^\infty f(x)dx = \int_0^1 f(a + t/(1-t)). 1 / (1-t)^2 dt
  ## as suggested on http://ab-initio.mit.edu/wiki/index.php/Cubature

  umax3 <- 1; umin3 <- 0;
  C3 <- (umax3 - umin3)/2 ; D3 <- (umax3 + umin3)/2
    
  unodes3 <- C3 * GL3$nodes + D3
  uweights3 <- GL3$weights * C3 * 1 / (1 - unodes3)^2
  unodes3 <- ucut + unodes3 / (1 - unodes3)
    
  Nu3 <- length(unodes3)
  
  in3 <- integrand2GL(u=unodes3,
                    d=d, lambda=lambda,
                    epsilon=epsilon, thickness=thickness)
  
  weights3 <- matrix(uweights3, nrow=Nlambda, ncol=Nu3, byrow=TRUE)
  
  integral3.perp <- rowSums(in3$integrand.p*weights3)
  integral3.par <- rowSums(in3$integrand.s*weights3)

  ## for Mrad, we use the same integration points as GL1 because we study the radiative region

  thetamax <- pi/2; thetamin <- 0;
  C4 <- (thetamax - thetamin)/2 ; D4 <- (thetamax + thetamin)/2
  thetanodes <- C4 * GL1$nodes + D4
  thetaweights <- GL1$weights * C4

  Ntheta <- length(thetanodes)
  
  in4 <- integrand.radGL(theta=thetanodes,
                      d=d, lambda=lambda,
                      epsilon=epsilon, thickness=thickness)
  
  weights4 <- matrix(thetaweights, nrow=Nlambda, ncol=Ntheta, byrow=TRUE)

  Mrad.perp <- 3/4 * rowSums(in4$integrand.p*weights4)
  Mrad.par <- 3/8 * rowSums(in4$integrand.s*weights4)
    
  data.frame(wavelength=lambda,
             Mtot.perp = 1 + 3/2*(integral1.perp + integral2.perp + integral3.perp),
             Mtot.par = 1 + 3/4*(integral1.par + integral2.par + integral3.par),
             Mrad.perp = Mrad.perp, Mrad.par = Mrad.par) 

}
