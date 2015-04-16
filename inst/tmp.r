# library(planar)
# 
# lambda <- 600.0
# theta <- pi/4
# epsilon <- c(incident = 1.5^2, 1.33^2)
# thickness <- c(0, 0)
# 
# a <- recursive_fresnel(wavelength = lambda, angle = theta, 
#                   epsilon = epsilon,
#                   thickness = thickness, polarisation = "p")
# 
# k0 <- 2*pi/lambda
# kx <- k0*sin(theta)
# planar$recursive_fresnel(as.vector(k0), as.matrix(kx), 
#                          as.matrix(epsilon), 
#                          as.vector(thickness), as.integer(0L))
# 
# planar$multilayer(as.vector(k0), as.matrix(kx), 
#                          as.matrix(t(epsilon)), 
#                          as.vector(thickness), as.integer(0L))

library(Rcpp)
library(RcppArmadillo)

# a <- cppFunction(depends = "RcppArmadillo",
#             '
# arma::mat euler(const double phi, const double theta, const double psi)
# {
#   arma::mat Rot(3,3);
#   const double cosphi = cos(phi), cospsi = cos(psi), costheta = cos(theta);
#   const double sinphi = sin(phi), sinpsi = sin(psi), sintheta = sin(theta);
#   Rot(0,0) = cospsi*cosphi - costheta*sinphi*sinpsi;
#   Rot(0,1) = cospsi*sinphi + costheta*cosphi*sinpsi; 
#   Rot(0,2) = sinpsi*sintheta;
#   
#   Rot(1,0) = -sinpsi*cosphi - costheta*sinphi*cospsi; 
#   Rot(1,1) = -sinpsi*sinphi + costheta*cosphi*cospsi; 
#   Rot(1,2) = cospsi*sintheta;
#   
#   Rot(2,0) = sinphi*sintheta;
#   Rot(2,1) = -cosphi*sintheta; 
#   Rot(2,2) = costheta;
#   return (Rot);
# }')
# 
# a(pi/2, 0, 0)

library(cda)
cda$euler(pi/2, 0.0, 0.0)

