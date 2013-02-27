## comparison of the calculation of near field enhancement outside of a thin metal film
## with Fresnel reflection and transmission coefficients

library(planar)
require(reshape2)

wvl <- seq(200, 1000,by=2)*1e-3
gold <- epsAu(wvl*1e3)

field.outside <- function(d=1, theta = seq(0,pi/2-0.001,length=500),
                          epsilon=list(incident = 1.0^2, gold$epsilon[1], 1.5^2),
                          thickness=c(0, 45, 0),
                          wavelength=633, polarisation="p", ...){

 res <- multilayer(wavelength, theta=theta, epsilon=epsilon, d=rep(d,length(thickness)),
                   thickness=thickness, polarisation=polarisation, ...)

 n <- length(epsilon)
 sinl <- sin(theta)
 sinr <- sqrt(epsilon[[1]]) / sqrt(epsilon[[n]]) * sinl

 k0 <- 2*pi/wavelength
 kx <- k0*sqrt(epsilon[[1]]) * sinl
 kzl <- sqrt(epsilon[[1]] * k0^2 - kx^2+0i)
 kzr <- sqrt(epsilon[[n]] * k0^2 - kx^2+0i)
 
 pol.fac1 <- if(polarisation == "p") epsilon[[1]] / epsilon[[n]]  else 1
 pol.fac2 <- if(polarisation == "p") sinl^2  else 1
 pol.fac3 <- if(polarisation == "p") sinr^2  else 1
 
 left <- pol.fac2 * Mod(exp(-1i*d* kzl) + res$reflection * exp(1i*d* kzl))^2
 right <- pol.fac3 * Mod(res$transmission * exp(1i*(d)* kzr))^2
 
 left2 <- if(polarisation == "p") res$Ml.perp[[1]][,1] else res$Ml.par[[1]][,1]
 right2 <- if(polarisation == "p") res$Mr.perp[[2]][,2]  else res$Mr.par[[2]][,2] 
 d <- data.frame(theta=theta*180/pi,
                 left=left, right=right*pol.fac1,
                 left2 = left2, right2 = right2)

 classify(d, id="theta", vars=list(side = rep(c("left", "right"), 2),
                           model = rep(c("matrix","fresnel"), each=2)))
 
}

test <- field.outside(10, thickness=c(0, 50, 0), polarisation="s",
                      epsilon=list(incident = 1.45^2, gold$epsilon[1], 1.0^2))

p <- 
ggplot(test) + 
  geom_path(aes(theta, value, colour=side, linetype=model), size=1.2)

p

