## comparison of the 4 codes to calculate R, T, A

library(planar)
require(reshape2)
require(ggplot2)

wvl <- seq(200, 1000,by=2)*1e-3
gold <- epsAu(wvl*1e3)

RTA_comparison <- function(angle = 40*pi/180,
                          epsilon=list(incident = 1.52^2, gold$epsilon, 1.1^2),
                          thickness=c(0, 45, 0), polarisation="p",
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
#  all <- lapply(all, transform, unity=1-R-T-A)
 melt(all, id="wavelength")
 
}

testp <- RTA_comparison()
tests <- RTA_comparison(pol="s")

p <- 
ggplot(testp, aes(wavelength, value, colour=L1, linetype=variable,group=interaction(L1, variable)))+
  geom_line(position=position_jitter(width=0.)) +
  scale_y_continuous( )

gridExtra::grid.arrange(p, p %+% tests)