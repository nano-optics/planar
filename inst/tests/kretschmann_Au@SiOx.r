## Reflectivity against internal incident angle
## for the Kretschmann configuration, at fixed wavelength

library(planar)
library(ggplot2)

wvl <- 632.8
silver <- epsAg(wvl)
gold <- epsAu(wvl)


model <- function(thickness=0){
  theta <- seq(30,60,length=2e3)
parameters <- list(epsilon=list(1.5^2, gold$epsilon, 1.5^2, 1.0),
                   lambda=gold$wavelength*1e3, thickness=c(0, 50, thickness, 0),
                   theta=theta*pi/180, polarisation='p')

d <- do.call(recursive.fresnel2, parameters)
data.frame(angle=theta, R = d$R)
}

params <- data.frame(thickness=seq(0, 20, by=5))
m <- mdply(params, model)

tir <- asin(sqrt(parametersAu$epsilon[[3]]) /
            sqrt(parametersAu$epsilon[[1]])) * 180/pi
p <- 
ggplot(m) +
  geom_vline(aes(xintercept=x),
             data=data.frame(x=tir),
             linetype=2,color="grey50") +
  geom_path(aes(angle, R, color=factor(thickness))) +
  labs(colour="thickness /nm") + scale_color_brewer(pal="Set1")+
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees") +
  theme_bw()

p 
ggsave("shift-SiOx.pdf",p)
ddply(m, .(thickness), summarize, spp = angle[which.min(R)])
## thickness      spp
## 1         0 44.42221
## 2         5 45.08254
## 3        10 45.83292
## 4        15 46.71836
## 5        20 47.72386
