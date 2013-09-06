## Reflectivity against internal incident angle
## for the Kretschmann configuration, at fixed wavelength

library(planar)
library(ggplot2)

wvl <- 632.8
silver <- epsAg(wvl)
gold <- epsAu(wvl)
gold$epsilon <- -11.40185 + 1i


model <- function(thickness=0){
  theta <- seq(40,50,length=2e3)
parameters <- list(epsilon=list(1.5^2, gold$epsilon, 1.5^2, 1.0),
                   lambda=gold$wavelength, thickness=c(0, 50, thickness, 0),
                   theta=theta*pi/180, polarisation='p')

d <- do.call(recursive_fresnel, parameters)
  print(str(d))
data.frame(angle=theta, R = d$R, q=d$q)
}

params <- data.frame(thickness=seq(0, 20, by=5))
m <- mdply(params, model)

p <- 
ggplot(m) +
  geom_line(aes(angle, R, colour=thickness, group=thickness)) +
  labs(colour="thickness /nm") + #scale_color_brewer(palette="Set1")+
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees") +
  theme_bw()

p 
# ggsave("shift-SiOx.pdf",p)
# ddply(m, .(thickness), summarize, spp = angle[which.min(R)])
## variable      spp
## 1     gold 44.43722
## 2   silver 43.58179
