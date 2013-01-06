## radiation pattern of a dipole near a dielectric/(metal)/dielectric interface
## parallel and perpendicular orientations, p- polarisation

library(planar)
require(reshape2)
library(ggplot2)
library(gridExtra)


## from left to right
## incident water | metal | glass

params <- function(wavelength, ntheta=300)
  list(lambda = wavelength,
       epsilon =list(1.5^2, epsAg(wavelength)$epsilon, 1.33^2),
                 theta = seq(70, 78, length=ntheta)[-c(1, ntheta)] * pi/180,
                thickness = c(0, 50, 0),
                d = 1,
                 polarisation = 'p')



raman_enhancement <- function(wavelength = 568, shift=0, ntheta=10){
  wvl <- raman.shift(wavelength, shift)
  p <- params(wvl, ntheta)
  Mstokes <- with(p,
            multilayer(lambda = lambda, theta = theta, polarisation = polarisation,  
                  thickness = thickness, d = d,
                  epsilon = epsilon))
  res <- Mstokes[['Mr.perp']][[2]]
  data.frame(theta = p[['theta']]*180/pi,
  Mperp = c(t(res)), shift = rep(shift, each=ncol(res)))
}

test <- raman_enhancement(shift = seq(-1000, 1000, length=200),ntheta=10)

p <- 
ggplot(test, aes(theta, Mperp, colour=shift)) +
  geom_line()

ggplot(test, aes(shift, Mperp, colour=theta, group=theta)) +
  geom_line() +
  guides(colour=guide_legend()) +
  labs(x = expression("Raman shift /"*cm^-1), colour="Angle") +
  theme_minimal()