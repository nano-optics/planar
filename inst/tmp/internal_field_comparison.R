library(planar)

lambda <- 633
k0 <- 2*pi/lambda
kx <- k0*sin(pi/4)*1.5
epsilon <- c(1.5^2, epsAg(lambda)$epsilon, 1.0^2, 1.0^2)
thickness <- c(0, 50, 10, 0)
psi <- 0
angle <- 0.7617526
angle <- 2*pi/180

struct <- tamm_stack(wavelength = 704,
                     lambda0=650, N=6, incidence = "left",
                     nH = 1.7, nL = 1.3, dm = 25, 
                     nleft = 1.0, nright = 1.52)


comp <- field_profile(wavelength = struct$wavelength, 
               angle = angle, 
               polarisation = "p", 
               thickness = struct$thickness, 
               dmax = 500, res = 100, 
               epsilon = struct$epsilon, displacement = FALSE)
              
d <- seq(min(comp$x), max(comp$x), length=300)
# d <- -10
k0 <- 2*pi/struct$wavelength
n1 <- sqrt(struct$epsilon[[1]])
res <- planar$multilayer_field(k0, k0*sin(angle)*n1, unlist(struct$epsilon),  
                               struct$thickness, d, psi)
field <- colSums(res$E * Conj(res$E)) 

plot(d, field, t="l")
require(plyr)
Mperp <- arrange(subset(comp, variable == "M.perp", select=c("x", "value")), x)
Mpar <- arrange(subset(comp, variable == "M.par", select=c("x", "value")), x)
lines(Mperp$x, Mpar$value +  Mperp$value, lty=2,col="red")

