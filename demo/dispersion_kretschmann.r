## Plot the dispersion curve of planar SPP on a thin metal film in the Kretschmann configuration
## the internal angle is varied artificially from 0 to 90 degrees

library(planar)
require(reshape2)


k0 <- seq(1e-4, 3e-2, length=500)
wvl <- 2*pi/k0
silver <- epsAg(wvl)
gold <- epsAu(wvl)

kretschmannAu <- recursive_fresnelcpp(k0=k0,
                                     q=seq(0,1, length=500),
                                     epsilon=list(1.5^2, gold$epsilon, 1.0),
                                     thickness=c(0, 50, 0),
                                     polarisation='p')

kretschmannAg <- recursive_fresnelcpp(k0=k0,
                                     q=seq(0,1, length=500),
                                     epsilon=list(1.5^2, silver$epsilon, 1.0),
                                     thickness=c(0, 50, 0),
                                     polarisation='p')
dispersion <- function(l, ...){
  
  m <- melt(data.frame(k0=l$k0, R=l$R), id=c("k0"))
  m$q <- rep(Re(l$q), each=nrow(l$R))
  invisible(m)
}

mAu <- dispersion(kretschmannAu)
mAu$material <- "Au"
mAg <- dispersion(kretschmannAg)
mAg$material <- "Ag"

m <- rbind(mAu,mAg)

# library(lattice)
# pal1 <- grey(seq(0,1,leng=100))
# p <- levelplot(value~q*k0|material, data=m, panel=panel.levelplot.raster,
#                interpolate=TRUE, layout=c(1,2),
#                scales=list(x=list(axs="i", at=seq(0,1,by=0.2)),
#                  y=list(relation="free")), xlim=c(0,1),
#                col.regions = colorRampPalette(pal1)(1e3), cut=1e3,
#                at=seq(0,1,length=1e3),
#                xlab = expression(q==k[x] / k[1]), ylab=expression(k[0]))
# 
# p

p <- 
ggplot(m, aes(q, k0, fill=value)) +
  facet_wrap(~material, ncol=1) +
  geom_raster() + labs(fill = "R") +
  scale_x_continuous(expression(q==k[x] / k[1]), expand=c(0,0))+
  scale_y_continuous(expression(k[0]/nm^-1), expand=c(0,0))+
  theme_minimal()

p
