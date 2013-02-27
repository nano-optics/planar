
## Reproducing Fig. 6.6, p. 188 of Mac Leod's Thin Film Optical Filters
## the structure is a stack of lambda/4 layers of indices nH and nL on a glass substrate
## with increasing number of layers, the reflectivity stop-band becomes stronger

library(planar)
library(ggplot2)
require(reshape2)
require(plyr)

makeStack <- function(n = 3, lambda=seq(200, 1000),
                      lambda0 = 460, thickness = lambda0/4, nH=2.3, nL=1.38, nS=1.52,
                      angle=0){

  epsilon.list <- c(1, rep(c(nL, nH)^2, n), nS)
  thickness.list <- c(0, rep(thickness/c(nL, nH), n), 0)
  
  params <- list(epsilon=as.list(epsilon.list),
                 lambda=lambda, thickness=thickness.list,
                 theta=angle*pi/180, polarisation='p')
  
  data.frame(do.call(recursive_fresnelcpp, params))
}

params <- expand.grid(n=seq(5,20,by=2), angle=seq(0,60, by=30))
all <- mdply(params, makeStack)

p <- 
ggplot(all)+ facet_grid(angle~.)+
  geom_line(aes(2*pi/k0, R, colour=n, group=n))+
  labs(colour="layers") +
  scale_x_continuous("wavelength /nm")+
  scale_y_continuous("Reflectivity", expand=c(0,0), lim=c(0,1))

p
