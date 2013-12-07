
library(planar)
library(ggplot2)
library(plyr)

results <- multilayer(epsilon=list(1.52^2, epsAg(632.8)$epsilon, 1.0),
                         wavelength=632.8, thickness=c(0, 50, 0), d=1,
                         angle=seq(0, pi/2, length=2e3), polarisation='p')

plot(results$angle, results$Mr.perp[[2]], t="l")
maxi <- max(results$Mr.perp[[2]], na.rm=T)
spp <- results$angle[which.min(results$R)]

simul <- function(w0=10){
  w0 <- w0*1e3
  xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0+5000,length=100), y=0, z=seq(0, 300, by=50)))
  res <- adply(xyz, 1, gaussian_field, d=50, nl=sqrt(epsAg(632.8)$epsilon),
               wavelength=632.8, w0=w0, ni=1.52, alpha=spp, maxEval=1000)
  data.frame(xyz, field=res[[2]])
}
# 

all <- mdply(data.frame(w0=c(1e3, 5e2, 1e2, 50, 10, 5, 1)), simul, .progress="text")

p <- 
  ggplot(all, aes(x/w0/1000, field, group=z, colour=z))+
  geom_line() + facet_grid(w0~., scales="free", labeller=label_both) +
  geom_vline(aes(x=0,y=NULL),lty=2) +
#   geom_hline(aes(x=0,yintercept=maxi),lty=3) +
  labs(x=expression(x/w[0]), y=expression("|E|"^2), 
       colour=expression("z /nm")) +
  coord_cartesian(xlim=c(-5,5)) + theme_minimal()+
  guides(colour=guide_legend())
print(p)

# library(animation)
# saveGIF({
#   ani.options(nmax = 30)
#   l_ply(rev(seq(1e4, 1e6, length=10)), simul, .progress="text")
# }, interval = 0.05, movie.name = "coupling.gif", ani.width = 600, ani.height = 600)
# 
# 
# 
