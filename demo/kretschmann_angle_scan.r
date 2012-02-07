## Reflectivity against internal incident angle
## for the Kretschmann configuration, at fixed wavelength

library(planar)
library(ggplot2)
require(reshape2)
require(plyr)

wvl <- 632.8
silver <- epsAg(wvl)
gold <- epsAu(wvl)

## with Cr layer
parametersAu <- list(epsilon=list(1.515^2, -1.23 +20.78i, gold$epsilon, 1.0^2),
                   lambda=gold$wavelength*1e3, thickness=c(0, 2, 50, 0),
                   theta=seq(0,pi/2,length=2e3), polarisation='p')

parametersAu <- list(epsilon=list(1.5^2, gold$epsilon, 1.0),
                   lambda=gold$wavelength*1e3, thickness=c(0, 50, 0),
                   theta=seq(0,pi/2,length=2e3), polarisation='p')


parametersAg <- list(epsilon=list(1.5^2, silver$epsilon, 1.0),
                   lambda=silver$wavelength*1e3, thickness=c(0, 50, 0),
                   theta=seq(0,pi/2,length=2e3), polarisation='p')

dAu <- do.call(recursive.fresnel2, parametersAu)
dAg <- do.call(recursive.fresnel2, parametersAg)

m <- melt(data.frame(angle = parametersAu$theta*180/pi,
                     gold = dAu$R, silver = dAg$R), id="angle")

tir <- asin(sqrt(parametersAu$epsilon[[length(parametersAu$epsilon)]]) /
            sqrt(parametersAu$epsilon[[1]])) * 180/pi
p <- 
ggplot(m) +
  geom_vline(aes(xintercept=x),
             data=data.frame(x=tir),
             linetype=2,color="grey50") +
  geom_path(aes(angle, value, color=variable)) +
  labs(colour="material") + scale_color_brewer(palette="Set1")+
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees", expand=c(0,0), limits=c(0,90),
                     breaks=sort(c(seq(0,90,by=15), round(tir,1)))) +
  theme_bw()

p 

ddply(m, .(variable), summarize, spp = angle[which.min(value)])
## variable      spp
## 1     gold 44.43722
## 2   silver 43.58179
