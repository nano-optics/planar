## Visualisation of the power loss as a function of in-plane wavevector for a dipole near a thin metal film
## note the divergence at the light line (q=1)
## and the slow decay of the quasi-static image dipole interaction at large q

library(planar)
library(ggplot2)
require(reshape2)
library(lattice)

wvl <- seq(200, 1000,by=2)
silver <- epsAg(wvl)
gold <- epsAu(wvl)

## focus on the radiative and SPP region
q <- sort(unique(c(seq(0.2,1, length=500), seq(1,1.2, length=500)))) 

integrand <- function(d=1, m= "silver", q, ...){
  mat <- get(m)
  int <- dipole.integrand(d=d, q=q, lambda= mat$wavelength*1e3,
                          epsilon = list(incident=1.0^2, mat$epsilon, 1.0^2),
                          thickness = c(0, 50, 0))

  ## scaling <- 10 # max(int$integrand.p) / max(int$integrand.s)
  m <- data.frame(wavelength = rep(mat$wavelength*1e3, length(q)),
                  k0=2*pi/rep(mat$wavelength*1e3, length(q)),
                  q = rep(q, each=length(mat$wavelength)),
                  perpendicular=c(int$integrand.p),
                  parallel=c(int$integrand.s))
  melt(m, id=1:3)
}

all <- integrand(q=q)


p <- levelplot(value~q*k0 | variable, data=all, panel=panel.levelplot.raster,
               interpolate=TRUE, colorkey=FALSE, layout=c(1,2),
               col.regions = colorRampPalette(grey(seq(0,1,length=10)))(1e3), cut=1e3,
               ## at= seq(0, max(m$value),length=1e3),
               scales=list(x=list(at=seq(0,1.2,by=0.2))),
               xlab = expression(q==k[x] / k[1]), ylab=expression(k[0]))

p

## now do slices for a wide view

q <- sort(unique(c(seq(0.5,0.9999, length=500), seq(1.001,1.2, length=500),
                   seq(1.2,1000, length=1000))))

wvl <- seq(0.2, 0.8, by=0.2)
silver <- epsAg(wvl*1e3)
all <- integrand(q=q, d=1)

## str(all)

p2 <- ggplot(all) + facet_grid(variable~.)+
  geom_path(aes(q, value, color=factor(wavelength))) +
  coord_cartesian(ylim=c(-60, 120))+
  scale_y_continuous(breaks=seq(-50, 120, by=50 ))+
  scale_x_log10() + labs(y="", colour=expression(lambda/nm))+
  theme_minimal()

p2
