##' Multilayer Fresnel coefficients
##'
##' computes the reflection coefficient of a multilayered interface
##' @title recursive.fresnel
##' @export
##' @param lambda [vector] wavelength in nm
##' @param k0 [vector] wavevector in nm^-1
##' @param theta [vector] incident angles in radians
##' @param q [vector] normalised incident in-plane wavevector
##' @param epsilon list of N+2 dielectric functions, each of length 1 or length(lambda)
##' @param thickness vector of N+2 layer thicknesses, first and last are dummy
##' @param polarisation [character] switch between p- and s- polarisation
##' @return fresnel coefficients and field profiles
##' @author baptiste Auguie
##' @examples
##' library(planar)
##' demo(package="planar")
recursive.fresnel <- function(lambda = NULL, k0 = 2*pi/lambda,
                       theta = NULL, q = sin(theta),
                       epsilon = list(incident=1.5^2, 1.33),
                       thickness = c(0, 0),
                       polarisation = c('p', 's')){

  ## checks
  stopifnot(thickness[1]==0L, thickness[length(thickness)]==0L)
  polarisation <- match.arg(polarisation)
  
  ## define constants
  Nlambda <- length(k0)
  Nq <- length(q)
  Nlayer <- length(thickness)
  k02 <- k0^2
  kx <- outer(k0*sqrt(epsilon[[1]] + 0i), q) # kx = q*k0
  kx2 <- kx^2
  
  ## loop to calculate kiz 
  kiz <- array(0 + 0i, dim=c(Nlambda, Nq, Nlayer))
  
  for (ii in seq(1, Nlayer)){
    kiz[ , , ii] <- sqrt(matrix(epsilon[[ii]]*k02, nrow=Nlambda, ncol=Nq) - kx2 + 0i)
  }

  ## calculate all single-interface fresnel coefficients and phase factors

  Nsurface <- Nlayer -1 
  rsingle <- tsingle <- array(0 + 0i, dim=c(Nlambda, Nq, Nsurface))
  phase1 <- phase2 <- array(1 + 0i, dim=c(Nlambda, Nq, Nlayer))
  
  for (ii in seq(1, Nsurface)){
    
   if(polarisation == 'p'){

     a <- kiz[,,ii] / epsilon[[ii]]
     b <- kiz[,,ii+1] / epsilon[[ii+1]]
     
   } else { # s-polarisation
     
     a <- kiz[,,ii] 
     b <- kiz[,,ii+1]

   }
   
   rsingle[,,ii] <-  (a - b) / (a + b)
   tsingle[,,ii] <-  2 * a / (a + b)

   phase1[,,ii] <- exp(1i*thickness[ii]*kiz[,,ii])
   phase2[,,ii] <- exp(2i*thickness[ii]*kiz[,,ii])
   
 }
  phase1[,,Nlayer] <- 1 + 0i # 0 thickness for last medium
  phase2[,,Nlayer] <- 1 + 0i # 0 thickness for last medium
  
  ## now recursion, r.tmp is the combined reflection
  ## {r12, 23, ...}, then {r13, r24, ...}
  ## finally {r1N}

  ## starting from rsingle
  r.tmp <- rsingle
  t.tmp <- tsingle

  
  for(ii in seq(2, Nsurface)){
    
    for(jj in seq(1, Nlayer - ii)){
      r.tmp[,,jj] <- ( rsingle[,,jj] + r.tmp[,,jj+1]*phase2[,,jj+1]) /
        (1 + rsingle[,,jj]*r.tmp[,,jj+1]*phase2[,,jj+1])
      t.tmp[,,jj] <- ( tsingle[,,jj] * t.tmp[,,jj+1]*phase1[,,jj+1]) /
        (1 + rsingle[,,jj]*r.tmp[,,jj+1]*phase2[,,jj+1])
    }
  }
  
 
  
  list(wavelength=lambda, k0=k0, q=q,
       reflection=r.tmp[,,1], transmission=t.tmp[,,1],
       R=Mod(r.tmp[,,1])^2, T=Mod(t.tmp[,,1])^2)
  
  
}

##' Multilayer Fresnel coefficients
##'
##' computes the reflection coefficient of a multilayered interface
##' @title recursive.fresnel2
##' @export
##' @param lambda [vector] wavelength in nm
##' @param k0 [vector] wavevector in nm^-1
##' @param theta [vector] incident angles in radians
##' @param q [vector] normalised incident in-plane wavevector
##' @param epsilon list of N+2 dielectric functions, each of length 1 or length(lambda)
##' @param thickness vector of N+2 layer thicknesses, first and last are dummy
##' @param polarisation [character] switch between p- and s- polarisation
##' @return fresnel coefficients and field profiles
##' @author baptiste Auguie
##' @examples
##' library(planar)
##' demo(package="planar")
recursive.fresnel2 <- function(lambda = NULL, k0 = 2*pi/lambda,
                       theta = NULL, q = sin(theta),
                       epsilon = list(incident=1.5^2, 1.33),
                       thickness = c(0, 0),
                       polarisation = c('p', 's')){
  
  kx <- outer(k0*sqrt(epsilon[[1]]), q) # kx = q*k0
  epsilon = do.call(cbind, epsilon)
  ## case pure scalars
  if(nrow(epsilon) == 1L)
    epsilon <- matrix(epsilon, nrow=length(k0), ncol=length(thickness), byrow=TRUE)
  
  polarisation = if(polarisation == "p") 0L else 1L

  ## checks
  stopifnot(thickness[1]==0L,
            thickness[length(thickness)]==0L)
  
  stopifnot(length(thickness) == ncol(epsilon),
            nrow(epsilon) == length(k0),
            nrow(kx) == length(k0),
            ncol(kx) == length(q))

  ## call the C++ function
  res <- planar$recursive_fresnel(as.vector(k0),
                                  as.matrix(kx),
                                  as.matrix(epsilon),
                                  as.vector(thickness),
                                  as.integer(polarisation))
  
  list(k0 = k0, theta=theta, q=q,
       reflection=drop(res$reflection), 
       transmission=drop(res$transmission), 
       R=Mod(drop(res$reflection))^2,
       T=Mod(drop(res$transmission))^2)
}
