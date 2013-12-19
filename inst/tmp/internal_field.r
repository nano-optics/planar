library(planar)


lambda <- 633
k0 <- 2*pi/lambda
kx <- k0*sin(pi/4)*1.5
epsilon <- c(1.5^2, epsAg(lambda)$epsilon, 1.1^2, 1.0^2)
thickness <- c(0, 50, 10, 0)
z <- 10
psi <- 0
res <- planar$multilayer_field(k0, kx, epsilon,  thickness, z=70, psi)


simul <- function(angle, psi=0, pol="p"){
  n1 <- sqrt(epsilon[[1]])
  res <- planar$multilayer_field(k0, k0*sin(angle)*n1, epsilon,  thickness, z, psi)
  if(pol == "p")
   return(Mod(res$rp)^2) else return(Mod(res$rs)^2)
}

theta <- seq(0, pi/2, length=300)
plot(theta, sapply(theta, simul), t="l")
lines(theta, sapply(theta, simul, pol="s"), lty=2)


simul <- function(angle, psi=0){
  n1 <- sqrt(epsilon[[1]])
  res <- planar$multilayer_field(k0, k0*sin(angle)*n1, epsilon,  thickness, 51, psi)
  field <- res$E
  crossprod(field, Conj(field))
}

theta <- seq(0, pi/2-0.00001, length=300)
plot(theta*180/pi, sapply(theta, simul), t="l")
lines(theta*180/pi, sapply(theta, simul, psi=pi/2), lty=2)


index <- function(d){
  if(d<0) return(1.5^2)
  if(d>0 && d <=50) return(epsAg(lambda)$epsilon)
  if(d>50 && d <=60) return(1.1^2)
  if(d>60) return(1.0^2)
  
}

simul <- function(d, psi=0){
  angle <- 0.7617526
  n1 <- sqrt(epsilon[1])
  res <- planar$multilayer_field(k0, k0*sin(angle)*n1, epsilon,  thickness, d, psi)
  field <- res$E
#   field[c(1,2)] <- 0
#   field[3] <- 0
  crossprod(field, Conj(field)) #* index(d)^2
}

td <- seq(-100, 500, length=300)
plot(td, sapply(td, simul), t="l")
lines(td, sapply(td, simul), lty=2)

# lines(theta*180/pi, sapply(theta, simul, psi=pi/2), lty=2)




# 
# simul <- function(){
# 
#   res <- planar$multilayer_field(k0, kx, epsilon,	thickness, z, psi)
# 
# 
#   
#   ## results
#   data.frame(wavelength=wavelength, k0 = k0, 
#              angle=angle, q=q, 
#              rp=rp, tp=tp, rs=rs, ts=ts,
#              Rp=Rp, Tp=Tp, Ap = 1 - Rp - Tp,
#              Rs=Rs, Ts=Ts, As = 1 - Rs - Ts,
#              Mperp=Mperp, Mpara=Mpara)
# }
# 
# # multilayer2(500, angle=0)
# 
# loop_lambda <- function(wavelength, ...){
#   
#   metal <- epsAu(wavelength)
#   multilayer2(wavelength, epsilon=c(1.5^2, metal$epsilon, 1.0^2), ...)$Rp
#   
# }
# 
# loop_angle <- function(angle, ...){
#   
#   metal <- epsAu(632.8)
#   multilayer2(632.8, epsilon=c(1.5^2, metal$epsilon, 1.33^2), angle=angle, ...)$Mpara
#   
# }
# ang <- seq(0, pi/2-0.001, length=300)
# # ang <- c(44,44.5)*pi/180
# comp <- multilayer(wavelength=632.8, angle=ang, d=1,
#                           epsilon=list(1.5^2, epsAu(632.8)$epsilon, 1.33^2),
#                           thickness=c(0,50,0), pol="p")
# coms <- multilayer(wavelength=632.8, angle=ang, d=1,
#                    epsilon=list(1.5^2, epsAu(632.8)$epsilon, 1.33^2),
#                    thickness=c(0,50,0), pol="s")
# # test <- sapply(seq(400, 800), loop_lambda, angle=44.5*pi/180)
# test <- sapply(ang, loop_angle,d=51)
# plot(ang,test, t="l")
# matlines(ang, comp$Mr.par[[2]]+coms$Mr.par[[2]], lty=2)
# 
