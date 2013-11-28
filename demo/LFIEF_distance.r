## comparison of the calculation of near field enhancement outside of a thin metal film
## with Fresnel reflection and transmission coefficients

library(planar)
require(reshape2)

wvl <- 600
gold <- epsAu(wvl)

enhancement_factor <- function(d=1, angle = seq(0,pi/2-0.001,length=500),
                          epsilon=list(incident = 1.0^2, 1.2 + 5i, gold$epsilon[1], 1.0^2),
                          thickness=c(0, 5, 45, 0),
                          wavelength=633, ...){

 res.p <- multilayer(wavelength, angle=angle, epsilon=epsilon, d=rep(d,length(thickness)),
                   thickness=thickness, polarisation="p", ...)
 res.s <- multilayer(wavelength, angle=angle, epsilon=epsilon, d=rep(d,length(thickness)),
                   thickness=thickness, polarisation="s", ...)
 
 res2.p <- recursive_fresnel(wavelength, angle=angle, epsilon=epsilon, 
                   thickness=thickness, polarisation="p", ...)
 res2.s <- recursive_fresnel(wavelength, angle=angle, epsilon=epsilon, 
                           thickness=thickness, polarisation="s", ...)

 n <- length(epsilon)
 sinl <- sin(angle)
 sinr <- sinl * sqrt(epsilon[[1]] / epsilon[[n]])
 
 k0 <- 2*pi/wavelength
 kx <- k0*sqrt(epsilon[[1]]) * sinl
 kzl <- sqrt(epsilon[[1]] * k0^2 - kx^2+0i)
 kzr <- sqrt(epsilon[[n]] * k0^2 - kx^2+0i)
 
 left.s <- Mod(exp(-1i*d* kzl) + res2.s$reflection * exp(1i*d* kzl))^2
 right.s <- Mod(res2.s$transmission * exp(1i*(d)* kzr))^2
 left.p <- sinl^2 * Mod(exp(-1i*d* kzl) + res2.p$reflection * exp(1i*d* kzl))^2
 right.p <- epsilon[[1]] / epsilon[[n]] *sinr^2 * Mod(res2.p$transmission * exp(1i*(d)* kzr))^2
 
 left2.p <- res.p$Ml.perp[[1]][,1] 
 left2.s <- res.s$Ml.par[[1]][,1]
 right2.p <- res.p$Mr.perp[[n-1]][,2]  
 right2.s <- res.s$Mr.par[[n-1]][,2] 
 
 d <- data.frame(angle=angle*180/pi,
                 left.p=left.p, right.p=right.p,
                 left.s=left.s, right.s=right.s,
                 left2.p=left2.p, right2.p=right2.p,
                 left2.s=left2.s, right2.s=right2.s)

 classify(d, id="angle", vars=list(side = rep(c("left", "right"), 4),
                                   polarisation = rep(c("p", "s"), each=2, length.out=8),
                                   model = rep(c("matrix","fresnel"), each=4)))
 
}

test <- enhancement_factor(1, thickness=c(0, 5, 50, 0), 
                      epsilon=list(incident = 1.45^2, -2 + 20i, gold$epsilon[1], 1.3^2))

p <- 
ggplot(test) + facet_wrap(~polarisation, scales="free") +
  geom_path(aes(angle, value, colour=side, linetype=model), size=1.2)

print(p)

