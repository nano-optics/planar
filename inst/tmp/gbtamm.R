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


require(tamm)

m <- ff_simulation(wavelength=seq(300, 1000),  
                   lambda0=650, N=6, incidence = "left",
                   nH = 1.7, nL = 1.3, dm = 25, 
                   nleft = 1.0, nright = 1.52)


p1 <- ggplot(m, aes(wavelength, R)) +
  geom_line() +
  labs(x="wavelength /nm", y="reflectivity")+
  theme_minimal() +
  theme(panel.background=element_rect())

min <- subset(m, R==min(R))


m <- nf_simulation(energy = l2e(min$wavelength),
                   lambda0=650, N=6, incidence = "left",
                   nH = 1.7, nL = 1.3, dm = 25, 
                   nleft = 1.0, nright = 1.52, dmax = 1000)

m <- recode_dbr(m) # add custom names for layers
m2 <- subset(m, variable == "M.par")
m2$x <- m2$x - 100
limits <- ddply(m2, .(L1), summarize, abb=unique(abb), 
                material=unique(material),
                xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-Inf, 
                ymax=max(m$value))

p2 <- 
  ggplot(m2) + 
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, 
                group=L1), fill="white", data=limits) +
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, 
                fill=material, group=L1), data=limits) +
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous(expression("|E|"^2),expand=c(0,0)) +
  geom_hline(yintercept=0) +
  geom_line(aes(x, value)) +
  scale_colour_manual("", values=col_palette) +
  scale_fill_manual("", values=fill_palette) +
  #   scale_fill_brewer("", palette="Pastel1") +
  #   scale_colour_brewer("", palette="Set1") +
  guides(colour="none",fill="none")+
  theme_minimal() 

struct <- tamm_stack(wavelength = min$wavelength,
                     lambda0=650, N=6, incidence = "left",
                     nH = 1.7, nL = 1.3, dm = 25, 
                     nleft = 1.0, nright = 1.52)
# # 
# struct <- tamm_stack(wavelength = min$wavelength,
#                      lambda0=650, N=1, incidence = "left",
#                      nH = 1.0, nL = 1.0, dm = 0, 
#                      nleft = 1.0, nright = 1.0)

# struct <- list(wavelength=600, epsilon=c(1.1^2, 1.1^2, 1.1^2),
#                thickness=c(0, 20, 0))
w0 <- 1e6
xyz <- as.matrix(expand.grid(x=0, 
                             y=0,
                             z=seq(-struct$wavelength, sum(struct$thickness)+struct$wavelength, 
                                   length=300)))
res <- adply(xyz, 1, gaussian_near_field2, epsilon=unlist(struct$epsilon), 
             thickness=struct$thickness, wavelength = struct$wavelength,
             alpha=0.0, w0=w0, maxEval=500, .progress="text")

xyz <- data.frame(xyz, field=res[[2]])

ggplot(xyz, aes(z, y=field))+
  geom_line() +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /um"), fill=expression("|E|"^2), 
       y=expression("z /um")) 
