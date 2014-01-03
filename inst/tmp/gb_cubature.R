## test the speed gain of using cubature at c++ level

library(planar)
library(plyr)
library(ggplot2)


angle <- 0.00001
# angle <- 45*pi/180
wavelength <- 632.8
metal <- epsAg(wavelength)$epsilon
epsilon <- list(1.5^2, metal, 1.0)
thickness <- c(0, 100, 0)
d <- thickness[2]
ni <- sqrt(epsilon[[1]])
nl <- sqrt(epsilon[[2]])
no <- sqrt(epsilon[[3]])

## first, check the plane wave result
results <- multilayer(epsilon=epsilon,
                      wavelength=wavelength, thickness=thickness, d=1,
                      angle=seq(0, pi/2, length=2e3), polarisation='p')

maxi <- max(results$Mr.perp[[2]] + results$Mr.par[[2]], na.rm=T)
spp <- results$angle[which.max(results$Mr.perp[[2]] + results$Mr.par[[2]])]
print(spp)

simulation <- function(w0=10){
  w0 <- w0*1e3
  xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0,length=100), y=0, z=thickness[2]+1))
  
  res <- gaussian_near_field2(xyz, wavelength=wavelength,
                              epsilon=unlist(epsilon), thickness=thickness,
                              w0=w0, alpha=spp, maxEval=1000, tol=1e-2, progress=TRUE)
  
  data.frame(xyz, field=res)
}

simulation2 <- function(w0=10){
  w0 <- w0*1e3
  xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0,length=100), y=0, z=thickness[2]+1))
  
  res <- adply(xyz, 1, gaussian_near_field, nl=nl, no=no, ni=ni,d=d,
               wavelength=wavelength, w0=w0, alpha=spp, maxEval=1000, tol=1e-2)
  data.frame(xyz, field=res[[2]])
  
}

params <- data.frame(w0=c(10, 50, 100, 1000, 5000))
system.time(all <- mdply(params, simulation2))

subset(all, field == max(field))

p <- ggplot(all, aes(x/w0, field, group=w0, colour=factor(w0)))+
  geom_line()  +
  geom_vline(aes(x=0,y=NULL),lty=2) +
  geom_hline(aes(x=0,yintercept=maxi),lty=3) +
  annotate("text", label="plane-wave", y=maxi, x=-2.5, vjust=0, fontface="italic") +
  labs(x=expression(x/w[0]), y=expression("|E|"^2), 
       colour=expression(w[0]/mu*m)) +
#   coord_cartesian(xlim=c(-5,5)) + theme_minimal()+
  guides(colour=guide_legend(reverse=TRUE)) +
  theme(panel.background=element_rect(fill=NA)) +
  theme()

print(p)
