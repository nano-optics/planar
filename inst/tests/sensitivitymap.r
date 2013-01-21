## Reflectivity against internal incident angle
## for the Kretschmann configuration, at fixed wavelength

library(planar)
library(ggplot2)
library(plyr)
wvl <- 632.8
silver <- epsAg(wvl)
gold <- epsAu(wvl)


model <- function(dn = 0.001, thickness=100){
  theta <- seq(0,90,length=2e3)

  eps1 <- list(1.5^2, gold$epsilon, (1.4)^2, 1.0)
  eps2 <- list(1.5^2, gold$epsilon, (1.4+dn)^2, 1.0)
  p1 <- p2 <- list(epsilon=eps1,
                   lambda=gold$wavelength, thickness=c(0, 50, thickness, 0),
                   theta=theta*pi/180, polarisation='p')
  p2[["epsilon"]] <- eps2
  d1 <- do.call(recursive_fresnelcpp, p1)
  d2 <- do.call(recursive_fresnelcpp, p2)

  data.frame(angle=theta, R = d1$R, dR = (d1$R - d2$R) / dn)
}

params <- expand.grid(thickness = c(10, 50, 100, 200, 300),
                     dn=seq(0.0001, 0.2, length=200))
m <- mdply(params, model, .progress="text")
library(RColorBrewer)
col <- brewer.pal(3,"PRGn")

p <- 
ggplot(m) + facet_wrap(~thickness, scales="free",ncol=1) +
  geom_raster(aes(angle, dn, fill=dR)) +
  labs(colour="thickness /nm") + 
  scale_fill_gradient2(expression(dR/dn),midpoint = 0,
                       low = col[1], mid=col[2], high=col[3])+
  scale_y_continuous("dn", expand=c(0,0))+
  scale_x_continuous("Internal angle /degrees", expand=c(0,0)) +
  theme_minimal()

ggsave("sensitivitymap.pdf", p, width=10, height=30)

shift <- ddply(m, .(thickness, dn), summarize, spp = angle[which.min(R)])

p <- qplot(dn,spp, colour=thickness, group=thickness, data=shift,geom="line")

