## Plot the dispersion curve of coupled-SPPs on a thin free standing metal film
## the internal angle is varied from 0 to 90 degrees
library(planar)
require(reshape2)
library(lattice)

wvl <- seq(200, 1000,by=2)
silver <- epsAg(wvl)
gold <- epsAu(wvl)

symmetricAu <- recursive.fresnel2(lambda=gold$wavelength*1e3,
                                  q = seq(1, 1.4, length=500),
                                  epsilon=list(1.5^2, gold$epsilon, 1.5^2),
                                  thickness=c(0, 50, 0),
                                  polarisation='p')

symmetricAg <- recursive.fresnel2(lambda=silver$wavelength*1e3,
                                  q = seq(1, 1.4, length=500),
                                  epsilon=list(1.5^2, silver$epsilon, 1.5^2),
                                  thickness=c(0, 50, 0),
                                  polarisation='p')

dispersion <- function(l, ...){
  
  m <- melt(data.frame(k0=l$k0, R=l$R), id=c("k0"))
  m$q <- rep(Re(l$q), each=nrow(l$R))
  invisible(m)
}

mAu <- dispersion(symmetricAu)
mAu$material <- "Au"
mAg <- dispersion(symmetricAg)
mAg$material <- "Ag"

m <- rbind(mAu,mAg)

pal1 <- grey(seq(0,1,leng=100))

p <- levelplot(value~q*k0|material, data=m, panel=panel.levelplot.raster,
               interpolate=TRUE, layout=c(1,2),
               scales=list(x=list(axs="i", at=seq(1,1.4,by=0.2)),
                 y=list(relation="free")), xlim=c(1,1.4),
               col.regions = colorRampPalette(pal1)(1e3), cut=1e3,
               at=seq(0,max(m$value),length=1e3),
               xlab = expression(q==k[x] / k[1]), ylab=expression(k[0]))

p

