
## ----, echo=FALSE,results='hide'-----------------------------------------
library(knitr)
library(ggplot2)
opts_chunk$set(fig.path="rtaconsistency/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=TRUE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----, results='hide'----------------------------------------------------
library(planar)
library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape2)

wvl <- seq(200, 1000,by=2)
gold <- epsAu(wvl)



## ----simulation----------------------------------------------------------
RTA_comparison <- function(angle = 30*pi/180,
                          epsilon=list(incident = 1.5^2, gold$epsilon,
                                       Re(gold$epsilon*0+ 1.0^2)),
                          thickness=c(0, 20, 0), polarisation="p",
                          wavelength=gold$wavelength, ...){

 res1 <- multilayer(wavelength, angle=angle, epsilon=epsilon, 
                   thickness=thickness, polarisation=polarisation, ...)
 res2 <- multilayercpp(wavelength, angle=angle, epsilon=epsilon, 
                    thickness=thickness, polarisation=polarisation, ...)
 res3 <- recursive_fresnel(wavelength, angle=angle, epsilon=epsilon, 
                   thickness=thickness, polarisation=polarisation, ...)
 res4 <- recursive_fresnelcpp(wavelength, angle=angle, epsilon=epsilon, 
                           thickness=thickness, polarisation=polarisation, ...)
 
 all <- list("multilayer" = data.frame(res1[c("wavelength","R", "T", "A")]),
             "multilayercpp" = data.frame(res2[c("wavelength","R", "T", "A")]),
             "fresnel" = data.frame(res3[c("wavelength","R", "T", "A")]),
             "fresnelcpp" = data.frame(res4[c("wavelength","R", "T", "A")]))
 m <-  melt(all, id="wavelength")
m 
}

testp <- RTA_comparison(polarisation="p")
tests <- RTA_comparison(polarisation="s")

p <- 
ggplot(testp, aes(wavelength, value, colour=L1, linetype=variable,
                  group=interaction(L1, variable)))+
  geom_line(position=position_jitter(width=0, height=0.)) +
  scale_y_continuous( )

grid.arrange(p, p %+% tests)


