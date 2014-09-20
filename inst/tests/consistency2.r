## comparison of the 4 codes to calculate R, T, A

library(planar)
require(reshape2)
require(ggplot2)

wvl <- seq(200, 1000,by=2)*1e-3
gold <- epsAu(wvl*1e3)


RTA_comparison <- function(stack, angle = 30*pi/180, polarisation="p",
                          wavelength=gold$wavelength, ...){

  epsilon <- stack$epsilon
  thickness <- stack$thickness
 res1 <- multilayer(wavelength, angle=angle, epsilon=epsilon, 
                   thickness=thickness, polarisation=polarisation, ...)
 res2 <- multilayercpp(wavelength, angle=angle, epsilon=epsilon, 
                    thickness=thickness, polarisation=polarisation, ...)
#  res2 <- res1
 res3 <- recursive_fresnel(wavelength, angle=angle, epsilon=epsilon, 
                   thickness=thickness, polarisation=polarisation, ...)
 res4 <- recursive_fresnelcpp(wavelength, angle=angle, epsilon=epsilon, 
                           thickness=thickness, polarisation=polarisation, ...)
#  res4 <- res3
 all <- list("multilayer" = data.frame(res1[c("wavelength","R", "T", "A")]),
             "multilayercpp" = data.frame(res2[c("wavelength","R", "T", "A")]),
             "fresnel" = data.frame(res3[c("wavelength","R", "T", "A")]),
             "fresnelcpp" = data.frame(res4[c("wavelength","R", "T", "A")]))
#  all <- lapply(all, transform, unity=1-R-T-A)
 m <-  melt(all, id="wavelength")
#  subset(m, !L1 %in% c("multilayercpp"))
#   subset(m, L1 %in% c("multilayer",  "multilayercpp"))
m 
}


stack <- list(thickness = c(0, 800, 0),
              epsilon = list(1.0^2, 1.5^2, 1.33^2))
stack <- list(thickness = c(0, 800, 0),
              epsilon = list(1.0^2, 1.5^2+0.1i, 1.33^2))
stack <- list(thickness = c(0, 800, 0),
              epsilon = list(1.0^2+0.1i, 1.5^2+0.1i, 1.33^2+0.1i))
test <- RTA_comparison(stack)

ggplot(test, aes(wavelength, value, colour=L1, linetype=variable,group=interaction(L1, variable)))+
  geom_line(position=position_jitter(width=0, height=0.)) +
  scale_y_continuous( )
