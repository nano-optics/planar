library(planar)
library(plyr)
library(ggplot2)

lambda <- 633
k0 <- 2*pi/lambda
kx <- k0*sin(pi/4)*1.5
epsilon <- c(1.5^2, epsAg(lambda)$epsilon, 1.0^2, 1.0^2)
thickness <- c(0, 50, 10, 0)
z <- 10
psi <- 0

gaussian$integrand_gb2(c(0.001, 2*pi/180), c(0,0,51), k0, psi, 10*pi/180, 5000, 
                     epsilon, thickness)



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


wavelength <- 632.8
metal <- epsAg(wavelength)$epsilon
## first, check the plane wave result
results <- multilayer(epsilon=list(1.5^2, metal, 1.0),
                      wavelength=wavelength, thickness=c(0, 50, 0), d=1,
                      angle=seq(0, pi/2, length=2e3), polarisation='p')

maxi <- max(results$Mr.perp[[2]] + results$Mr.par[[2]], na.rm=T)
spp <- results$angle[which.max(results$Mr.perp[[2]] + results$Mr.par[[2]])]

simulation <- function(w0=10){
  w0 <- w0*1e3
  xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0+5000,length=100), y=0, z=c(51)))
  res <- adply(xyz, 1, gaussian_near_field2, w0=w0, alpha=0.75986, maxEval=500)
  data.frame(xyz, field=res[[2]])
}

params <- data.frame(w0=c(2, 5, 10, 25, 50, 1e2, 500, 1e3, 1e4 ))
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



lambda <- 633
epsilon <- c(1.5^2, epsAu(lambda)$epsilon, 1.0^2, 1.0^2)
thickness <- c(0, 50, 10, 0)


w0 <- 5000
xyz <- as.matrix(expand.grid(x=seq(-15e3, 15e3, length=100), 
                             y=0,
                             z=seq(-5000, 500, length=100)))
res <- adply(xyz, 1, gaussian_near_field2, alpha=0.7732184, w0=w0, maxEval=500, .progress="text")
xyz <- data.frame(xyz, field=res[[2]])

ggplot(xyz, aes(x/1e3, z/1e3, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /um"), fill=expression("|E|"^2), 
       y=expression("z /um")) +
  coord_fixed()



lambda <- 633
epsilon <- c(1.5^2, epsAu(lambda)$epsilon, 1.1^2, 1.0^2)
thickness <- c(0, 50, 10, 0)
w0 <- 5000
xyz <- as.matrix(expand.grid(x=seq(-3e3, 15e3, length=100), 
                             y=0,
                             z=seq(-50, 100, length=100)))
res <- adply(xyz, 1, gaussian_near_field2, epsilon=epsilon, thickness=thickness,
             alpha=0.7732184, w0=w0, maxEval=500, .progress="text")
xyz <- data.frame(xyz, field=res[[2]])

ggplot(xyz, aes(x/1e3, z/1e3, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /um"), fill=expression("|E|"^2), 
       y=expression("z /um")) 

library(tamm)

tamm_stack(wavelength = e2l(energy), angle = 0,
           polarisation = "p", incidence = "left", w0 = 1.3,
           lambda0 = e2l(w0), N = 15, dH = lambda0/4,
           dL = lambda0/4, nH = 3.7, nL = 3, dm = 100, dxL = 0,
           dxH = 0, metal = epsAu(wavelength), nleft = nH,
           nright = 1, nbleft = nleft, nbright = nright,
           dbleft = 100, dbright = 100, ...)


lambda <- 633
epsilon <- c(1.5^2, 1.0^2, 1.0^2, 1.0^2)
thickness <- c(0, 100, 300, 0)
w0 <- 1000
xyz <- as.matrix(expand.grid(x=seq(-25e3, 25e3, length=100), 
                             y=0,
                             z=seq(-10e3, 500, length=100)))
res <- adply(xyz, 1, gaussian_near_field2, epsilon=epsilon, thickness=thickness,
             alpha=0*0.7732184+80*pi/180, w0=w0, maxEval=500, .progress="text")
xyz <- data.frame(xyz, field=res[[2]])

ggplot(xyz, aes(x/1e3, z/1e3, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /um"), fill=expression("|E|"^2), 
       y=expression("z /um")) 
