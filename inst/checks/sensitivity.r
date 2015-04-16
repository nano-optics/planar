## Reflectivity against internal incident angle
## for the Kretschmann configuration, at fixed wavelength

library(planar)
library(ggplot2)
library(plyr)
wvl <- 632.8
silver <- epsAg(wvl)
gold <- epsAu(wvl)


model <- function(dn = 0, thickness=100){
  theta <- seq(30,90,length=2e3)
parameters <- list(epsilon=list(1.5^2, gold$epsilon, (1.4+dn)^2, 1.0),
                   lambda=gold$wavelength, thickness=c(0, 50, thickness, 0),
                   theta=theta*pi/180, polarisation='p')

d <- do.call(recursive_fresnelcpp, parameters)
data.frame(angle=theta, R = d$R)
}

params <- data.frame(dn=seq(0, 0.1, length=10))
m <- mdply(params, model)

ggplot(m) +
  geom_path(aes(angle, R, color=factor(dn))) +
  labs(colour="thickness /nm") +# scale_color_brewer(pal="Set1")+
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees") +
  theme_bw()

shift <- ddply(m, .(dn), summarize, spp = angle[which.min(R)])

qplot(dn,spp,data=shift,geom="line")

