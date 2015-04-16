## Reflectivity against internal incident angle
## for the Kretschmann configuration, at fixed wavelength
## n = 3.138, k = 3.316.

library(planar)
library(ggplot2)

wvl <- 632.8
silver <- epsAg(wvl)
gold <- epsAu(wvl)
Cr <- (3.138 + 3.316i)^2
Ti <-  (2.705 + 3.767i)^2

model <- function(thickness=0){
  theta <- seq(30,60,length=2e3)
parameters <- list(epsilon=list(1.5^2, Cr, gold$epsilon, 1.0),
                   lambda=gold$wavelength*1e3, thickness=c(0,thickness, 50,  0),
                   theta=theta*pi/180, polarisation='p')

d <- do.call(recursive.fresnel2, parameters)
data.frame(angle=theta, R = d$R)
}

params <- data.frame(thickness=seq(0, 10, by=2))
m <- mdply(params, model)

## tir <- asin(sqrt(parametersAu$epsilon[[4]]) /
##             sqrt(parametersAu$epsilon[[1]])) * 180/pi
p <- 
ggplot(m) +
  ## geom_vline(aes(xintercept=x),
  ##            data=data.frame(x=tir),
  ##            linetype=2,color="grey50") +
  geom_path(aes(angle, R, color=factor(thickness))) +
  labs(colour="thickness /nm") + scale_color_brewer(pal="Set1")+
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees") +
  theme_bw()

p 
## ggsave("shift-Cr.pdf",p)
## ggsave("shift-Ti.pdf",p)
ddply(m, .(thickness), summarize, spp = angle[which.min(R)])

##  thickness      spp
## 1         0 44.42221
## 2         2 44.43722
## 3         4 44.43722
## 4         6 44.43722
## 5         8 44.45223
## 6        10 44.45223
