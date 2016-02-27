library(planar)
library(ggplot2)
library(plyr)

theme_set(theme_minimal() + theme(panel.background=element_rect(fill=NA)))

# 
# r2 <- c(0,0,51)
# lambda <- 633
# k0 <- 2*pi/lambda
# psi <- 0
# epsilon <- c(1.5^2, epsAu(lambda)$epsilon,1^2)
# thickness <- c(0, 50, 0)
# 
# foo <- function(angle=0){
#   rt <- c(angle*pi/180,0)
#   planar:::collection$integrand_collection(rt, r2, k0, psi, epsilon, thickness)
# }
# 
# theta <- seq(0,90, length=1e3)
# 
# exact <- function(angle){
# res <- cpp_multilayer_field(k0, 
#                                  k0 * sin(angle*pi/180) * sqrt(epsilon[1]),
#                                  epsilon, 
#                                  thickness, 51, psi)$E
# crossprod(res, Conj(res))
# }
# 
# plot(theta, sapply(theta, foo), t="l")
# lines(theta, sapply(theta, exact), col="red", lty=2)

## line example
xyz <- as.matrix(expand.grid(x=0, 
                             y=0, z=seq(0, 200, length=100)))
test <- collection_ml(xyz, psi=0, epsilon= c(1.5^2, epsAu(633)$epsilon,1^2), 
                      thickness=c(0,50,0),
                      omega=c(4,45)*pi/180, wavelength=633,
                      maxEval=2000.0)

xyz2 <- data.frame(xyz, I=drop(test))
ggplot(xyz2, aes(z, I, colour=x, group=x))+
  geom_line()


xyz <- as.matrix(expand.grid(x=2*seq(-200,200,length=30), 
                             y=2*seq(-200,200,length=30), z=1))

test <- collection_ml(xyz, omega=c(0,50)*pi/180, maxEval=500.0)


xyz2 <- data.frame(xyz, I=drop(test))
ggplot(xyz2, aes(x, y, fill= I))+
  geom_raster()


xyz <- as.matrix(expand.grid(x=2*seq(-200,200,length=30), 
                             z=2*seq(-500,300,length=30), y=0))

test <- collection_ml(xyz, omega=c(0,80)*pi/180, 
                      maxEval=500.0, psi=0)


xyz <- data.frame(xyz, I=drop(test))
ggplot(xyz, aes(x, z, fill= I))+
  geom_raster()


