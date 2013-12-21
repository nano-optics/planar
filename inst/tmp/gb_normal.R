library(planar)
library(plyr)
library(ggplot2)

gaussian_near_field2 <- function(x=1, y=1, z=1, wavelength=632.8, alpha = 15*pi/180, psi=0, 
                                 w0=1e4, epsilon = c(1.5^2, epsAg(lambda)$epsilon, 1.0^2, 1.0^2),
                                 thickness = c(0, 50, 10, 0),
                                 cutoff = min(1, 3*wavelength/(Re(sqrt(epsilon[1]))*pi*w0)), 
                                 maxEval = 3000, tol=1e-04, field=FALSE){
  
  k0 <- 2*pi/wavelength
  res <- cubature::adaptIntegrate(gaussian$integrand_gb2,
                                  lowerLimit=c(0, 0), # rho in [0,1], angle in [0,2*pi]
                                  upperLimit=c(cutoff, 2*pi), 
                                  fDim = 6, tol = tol,
                                  maxEval = maxEval,
                                  r2 = c(x, y, z), k0=k0, psi= psi, alpha=alpha,
                                  w0=w0, epsilon=epsilon, thickness=thickness)$integral
  
  E <- complex(real = res[1:3], imaginary=res[4:6])
  if(field) return(E)
  
  Re(crossprod(E, Conj(E)))
}


angle <- 0 #45*pi/180
wavelength <- 632.8
metal <- epsAg(wavelength)$epsilon
epsilon <- list(1.5^2, metal, 1.0)
thickness <- c(0, 50, 0)
## first, check the plane wave result
results <- multilayer(epsilon=epsilon,
                      wavelength=wavelength, thickness=thickness, d=1,
                      angle=angle+0*seq(0, pi/2, length=2e3), polarisation='p')

maxi <- max(results$Mr.perp[[2]] + results$Mr.par[[2]], na.rm=T)
spp <- results$angle[which.max(results$Mr.perp[[2]] + results$Mr.par[[2]])]

simulation <- function(w0=10){
  w0 <- w0*1e3
  xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0+5000,length=100), y=0, z=c(51)))
  res <- adply(xyz, 1, gaussian_near_field2, tol=1e-6, #cutoff=2000/w0,
               epsilon=unlist(epsilon), thickness=thickness, wavelength=wavelength,
               w0=w0, alpha=angle, maxEval=1000)
  data.frame(xyz, field=res[[2]])
}

params <- data.frame(w0=c(1, 10, 100, 1e3))
all <- mdply(params, simulation, .progress="text")

p <- ggplot(all, aes(x/w0/1000, field, group=w0, colour=factor(w0)))+
  geom_line()  +
  geom_vline(aes(x=0,y=NULL),lty=2) +
  geom_hline(aes(x=0,yintercept=maxi),lty=3) +
  annotate("text", label="plane-wave", y=maxi, x=-2.5, vjust=1, fontface="italic") +
  labs(x=expression(x/w[0]), y=expression("|E|"^2), 
       colour=expression(w[0]/mu*m)) +
  coord_cartesian(xlim=c(-5,5)) + theme_minimal()+
  guides(colour=guide_legend(reverse=TRUE)) +
  theme(panel.background=element_rect(fill=NA)) +
  theme()

print(p)
