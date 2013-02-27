## Plot the dispersion curve of coupled-SPPs on a thin free standing metal film

library(planar)
require(reshape2)


k0 <- seq(1e-4, 2e-2, length=500)
wvl <- 2*pi/k0
silver <- epsAg(wvl)
gold <- epsAu(wvl)

kretschmannAu <- recursive_fresnelcpp(k0=k0,
                                      q=seq(1,1.4, length=500),
                                      epsilon=list(1.5^2, gold$epsilon, 1.5^2),
                                      thickness=c(0, 50, 0),
                                      polarisation='p')

kretschmannAg <- recursive_fresnelcpp(k0=k0,
                                      q=seq(1,1.4, length=500),
                                      epsilon=list(1.5^2, silver$epsilon, 1.5^2),
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

p <- 
  ggplot(m, aes(q, k0, fill=log10(value))) +
  facet_wrap(~material, ncol=1) +
  geom_raster() + labs(fill = "log(R)")+
  scale_x_continuous(expression(q==k[x] / k[1]), expand=c(0,0))+
  scale_y_continuous(expression(k[0]/nm^-1), expand=c(0,0))+
  theme_minimal()

p
