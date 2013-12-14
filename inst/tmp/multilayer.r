multilayer2 <- function(wavelength = 2*pi/k0, k0 = 2*pi/wavelength,
                       angle = asin(q), q = sin(angle),
                       epsilon = c(1.5^2, 1, 1.33^2),
                       thickness = c(0, 50, 0), d=1){
  
  ## checks
  Nlayer <- 3
  stopifnot(thickness[1]==0L, thickness[length(thickness)]==0L)
  
  k02 <- k0^2
  kn2 <- k02*epsilon
  kx <- q*sqrt(epsilon[1])*k0
  kx2 <- kx^2
  
  ## calculate kiz 
  kiz <- sqrt(kn2 - kx2 +0i)
  
  ## calculate the transition matrix M

  Ms11 <- Ms22 <- 1 + 0i
  Ms21 <- Ms12 <- 0 + 0i
  
  Mp11 <- Mp22 <- 1 + 0i
  Mp21 <- Mp12 <- 0 + 0i
  
  M11new <- M22new <- M21new <- M12new <- 0 + 0i
  Mip11 <- Mip12 <- Mip21 <- Mip22 <- rep(0 + 0i, Nlayer-1)
  Mis11 <- Mis12 <- Mis21 <- Mis22 <- rep(0 + 0i, Nlayer-1)
  
  for (ii in seq(1, Nlayer-1)){
    
    Ksi <- kiz[ii+1] / kiz[ii]
    Kpi <- epsilon[ii] / epsilon[ii+1] * Ksi
    
    phasei <- exp(1i*thickness[ii]*kiz[ii])
    
    Mip11[ii] <- 0.5*(1+Kpi) / phasei
    Mip21[ii] <- 0.5*(1-Kpi) * phasei
    Mip12[ii] <- 0.5*(1-Kpi) / phasei
    Mip22[ii] <- 0.5*(1+Kpi) * phasei
    M11new <- Mp11*Mip11[ii] + Mp12*Mip21[ii]
    M21new <- Mp21*Mip11[ii] + Mp22*Mip21[ii]
    M12new <- Mp11*Mip12[ii] + Mp12*Mip22[ii]
    M22new <- Mp21*Mip12[ii] + Mp22*Mip22[ii]
    Mp11 <- M11new
    Mp12 <- M12new
    Mp21 <- M21new
    Mp22 <- M22new
    
    Mis11[ii] <- 0.5*(1+Ksi) / phasei
    Mis21[ii] <- 0.5*(1-Ksi) * phasei
    Mis12[ii] <- 0.5*(1-Ksi) / phasei
    Mis22[ii] <- 0.5*(1+Ksi) * phasei
    M11new <- Ms11*Mis11[ii] + Ms12*Mis21[ii]
    M21new <- Ms21*Mis11[ii] + Ms22*Mis21[ii]
    M12new <- Ms11*Mis12[ii] + Ms12*Mis22[ii]
    M22new <- Ms21*Mis12[ii] + Ms22*Mis22[ii]
    Ms11 <- M11new
    Ms12 <- M12new
    Ms21 <- M21new
    Ms22 <- M22new
    
    
  }
  
  ## calculate the Fresnel coefficients
  tp <- 1 / Mp11
  ts <- 1 / Ms11
  rp <- Mp21 * tp
  rs <- Ms21 * ts
  
  ## calculate the fields
  Hiy.H1y <- Hpiy.H1y <- Eix.E1 <- Epix.E1 <- Eiz.E1 <- Epiz.E1 <- rep(0+0i, Nlayer)
  Eiy.E1y <- Epiy.E1y <- Hix.H1 <- Hpix.H1 <- Hiz.H1 <- Hpiz.H1 <- rep(0+0i, Nlayer)
  
  ## p-polarisation
  
  Hiy.H1y[Nlayer] <- tp
  Hpiy.H1y[Nlayer] <- 0i
  
  AuxE1 <- sqrt(epsilon[1] + 0i) / k0 / epsilon[Nlayer]
  Eix.E1[Nlayer] <- Hiy.H1y[Nlayer] * kiz[Nlayer] * AuxE1
  Epix.E1[Nlayer] <- 0i
  AuxE2 <- epsilon[1] / epsilon[Nlayer] * Re(q)
  Eiz.E1[Nlayer] <- - Hiy.H1y[Nlayer] * AuxE2
  Epiz.E1[Nlayer] <- 0i   
  
  ## s-polarisation
  
  Eiy.E1y[Nlayer] <- ts
  Epiy.E1y[Nlayer] <- 0i
  
  AuxH1 <- 1 / (k0 * sqrt(epsilon[1] + 0i))
  Hix.H1[Nlayer] <- - Eiy.E1y[Nlayer] * kiz[Nlayer] * AuxH1
  Hpix.H1[Nlayer] <- 0i
  AuxH2 <- Re(q)
  Hiz.H1[Nlayer] <- Eiy.E1y[Nlayer] * AuxH2
  Hpiz.H1[Nlayer] <- 0i
  
  ## loop downwards to compute all field amplitudes
  for (ii in seq(Nlayer-1, 1, by=-1)){
    
    ## p-pol
    Hiy.H1y[ii] <- Mip11[ii]*Hiy.H1y[ii+1] + Mip12[ii]*Hpiy.H1y[ii+1]
    Hpiy.H1y[ii] <- Mip21[ii]*Hiy.H1y[ii+1] + Mip22[ii]*Hpiy.H1y[ii+1]
    
    AuxE1 <- sqrt(epsilon[1] + 0i) / k0 / epsilon[ii]
    Eix.E1[ii] <- Hiy.H1y[ii] * kiz[ii] * AuxE1
    Epix.E1[ii] <- - Hpiy.H1y[ii] * kiz[ii] * AuxE1
    AuxE2 <- epsilon[1] / epsilon[ii] * Re(q)
    Eiz.E1[ii] <- - Hiy.H1y[ii] * AuxE2
    Epiz.E1[ii] <- - Hpiy.H1y[ii] * AuxE2
    
    ## s-pol
    Eiy.E1y[ii] <- Mis11[ii]*Eiy.E1y[ii+1] + Mis12[ii]*Epiy.E1y[ii+1]
    Epiy.E1y[ii] <- Mis21[ii]*Eiy.E1y[ii+1] + Mis22[ii]*Epiy.E1y[ii+1]
    
    Hix.H1[ii]  <- - Eiy.E1y[ii] * kiz[ii] * AuxH1
    Hpix.H1[ii] <-   Epiy.E1y[ii] * kiz[ii] * AuxH1
    Hiz.H1[ii]  <-   Eiy.E1y[ii] * AuxH2
    Hpiz.H1[ii] <-   Epiy.E1y[ii] * AuxH2
  }
  
  ## right of interface 1
  if(d < 0 ){
    
    Mperp <- Mod(Eiz.E1[1]  * exp( 1i*d*kiz[1]) +
                   Epiz.E1[1] * exp(-1i*d*kiz[1]))^2
    Mparap <- Mod(Eix.E1[1]  * exp( 1i*d*kiz[1]) +
                    Epix.E1[1] * exp(-1i*d*kiz[1]))^2
    Mparas <- Mod(Eiy.E1y[1]  * exp( 1i*d*kiz[1]) +
                    Epiy.E1y[1] * exp(-1i*d*kiz[1]))^2
    Mpara <- Mparas + Mparap
    
  } else  if(d > 0 && d < thickness[2]){
    
    Mperp <- Mod(Eiz.E1[2]  * exp( 1i*d*kiz[2]) +
                   Epiz.E1[2] * exp(-1i*d*kiz[2]))^2
    Mparap <- Mod(Eix.E1[2]  * exp( 1i*d*kiz[2]) +
                    Epix.E1[2] * exp(-1i*d*kiz[2]))^2
    Mparas <- Mod(Eiy.E1y[2]  * exp( 1i*d*kiz[2]) +
                    Epiy.E1y[2] * exp(-1i*d*kiz[2]))^2
    Mpara <- Mparas + Mparap
    
  } else if(d > thickness[2]){
    d <-  thickness[2] -  d
    Mperp <- Mod(Eiz.E1[3]  * exp( 1i*d*kiz[3]) +
                   Epiz.E1[3] * exp(-1i*d*kiz[3]))^2
    Mparap <- Mod(Eix.E1[3]  * exp( 1i*d*kiz[3]) +
                    Epix.E1[3] * exp(-1i*d*kiz[3]))^2
    Mparas <- Mod(Eiy.E1y[3]  * exp( 1i*d*kiz[3]) +
                    Epiy.E1y[3] * exp(-1i*d*kiz[3]))^2
    Mpara <- Mparas + Mparap
    
  } 
  # ratio of refractive indices
  index.ratio <- Re(sqrt(epsilon[1])/sqrt(epsilon[Nlayer]))
  # ratio of cosines
  m <- Re(sqrt(1 - (index.ratio * q)^2 + 0i)/sqrt(1 - q^2 + 0i))
  
  rhop <- index.ratio
  rhos <- 1 / index.ratio
  
  Rp <- Mod(rp)^2
  Tp <- rhop * m * Mod(tp)^2
  Rs <- Mod(rs)^2
  Ts <- rhos * m * Mod(ts)^2
  
  ## results
  data.frame(wavelength=wavelength, k0 = k0, 
             angle=angle, q=q, 
             rp=rp, tp=tp, rs=rs, ts=ts,
             Rp=Rp, Tp=Tp, Ap = 1 - Rp - Tp,
             Rs=Rs, Ts=Ts, As = 1 - Rs - Ts,
             Mperp=Mperp, Mpara=Mpara)
}

# multilayer2(500, angle=0)

loop_lambda <- function(wavelength, ...){
  
  metal <- epsAu(wavelength)
  multilayer2(wavelength, epsilon=c(1.5^2, metal$epsilon, 1.0^2), ...)$Rp
  
}

loop_angle <- function(angle, ...){
  
  metal <- epsAu(632.8)
  multilayer2(632.8, epsilon=c(1.5^2, metal$epsilon, 1.33^2), angle=angle, ...)$Mpara
  
}

library(planar)
ang <- seq(0, pi/2-0.001, length=300)
# ang <- c(44,44.5)*pi/180
comp <- multilayer(wavelength=632.8, angle=ang, d=1,
                          epsilon=list(1.5^2, epsAu(632.8)$epsilon, 1.33^2),
                          thickness=c(0,50,0), pol="p")
coms <- multilayer(wavelength=632.8, angle=ang, d=1,
                   epsilon=list(1.5^2, epsAu(632.8)$epsilon, 1.33^2),
                   thickness=c(0,50,0), pol="s")
# test <- sapply(seq(400, 800), loop_lambda, angle=44.5*pi/180)
test <- sapply(ang, loop_angle,d=51)
plot(ang,test, t="l")
matlines(ang, comp$Mr.par[[2]]+coms$Mr.par[[2]], lty=2)

