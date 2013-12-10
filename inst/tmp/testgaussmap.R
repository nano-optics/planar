
library(planar)
library(ggplot2)
library(plyr)

results <- multilayer(epsilon=list(1.5^2, epsAu(632.8)$epsilon, 1.0),
                         wavelength=632.8, thickness=c(0, 50, 0), d=1,
                         angle=seq(0, pi/2, length=2e3), polarisation='p')

plot(results$angle, results$Mr.perp[[2]]+ results$Mr.par[[2]], t="l")
maxi <- max(results$Mr.perp[[2]] + results$Mr.par[[2]], na.rm=T)
spp <- results$angle[which.max(results$Mr.perp[[2]] + results$Mr.par[[2]])]
# 
# simul <- function(w0=10){
#   w0 <- w0*1e3
#   xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0+5000,length=100), y=0, z=seq(1)))
#   res <- adply(xyz, 1, gaussian_field, d=50, nl=sqrt(epsAu(632.8)$epsilon),
#                wavelength=632.8, w0=w0, ni=1.0, alpha=0, maxEval=2000)
#   data.frame(xyz, field=res[[2]])
# }
# # 
simul <- function(w0=10){
  w0 <- w0*1e3
  xyz <- as.matrix(expand.grid(x=seq(-15000, 15000,length=100), y=0, z=seq(-5e3, 500, length=100)))
  res <- adply(xyz, 1, gaussian_field, d=50, nl=sqrt(epsAu(632.8)$epsilon),
               wavelength=632.8, w0=w0, ni=1.5, alpha=0.7732184, maxEval=1000, .progress="text")
  data.frame(xyz, field=res[[2]])
}
# 
params <- data.frame(w0=2)
all <- mdply(params, simul)

ggplot(all, aes(x, z, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /nm"), fill=expression("|E|"^2), 
       y=expression("z /nm"))