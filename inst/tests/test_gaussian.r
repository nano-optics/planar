library(planar)
library(ggplot2)
library(plyr)

theme_set(theme_minimal() + theme(panel.background=element_rect(fill=NA)))

## line example
xyz <- as.matrix(expand.grid(x=seq(5e3,5e4,length=10), y=0, z=seq(1, 500, length=50)))
res <- adply(xyz, 1, gaussian_field, lambda=500, w0=1e4, ni=1.5+0i, alpha=60*pi/180,  .progress="text")
xyz <- data.frame(xyz, field=res[[2]])
ggplot(xyz, aes(z, field, colour=x, group=x))+
  geom_line()


w0 <- 1e3
## line example
xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0,length=200), y=0, z=seq(0, 100,by=10)))
res <- adply(xyz, 1, gaussian_field, lambda=500, psi=pi/2, maxEval = 1000, 
             w0=w0, ni=1.5+0i, alpha=70*pi/180,  .progress="text")
xyz <- data.frame(xyz, field=res[[2]])

ggplot(xyz, aes(x, field, colour=z, group=z))+
  geom_line() +
  geom_vline(aes(x=0,y=NULL),lty=2)

w0 <- 1e3
## line example
xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0,length=200), y=0, z=seq(0, 100,by=10)))
res <- adply(xyz, 1, gaussian_field, lambda=500, psi=pi/2, maxEval = 1000, 
             w0=w0, ni=1.5+0i, alpha=50*pi/180,  .progress="text")
xyz <- data.frame(xyz, field=res[[2]])
ggplot(xyz, aes(x, field, colour=z, group=z))+
  geom_line()


w0 <- 1e3
## 2D map
xyz <- as.matrix(expand.grid(x=seq(-3*w0, 3*w0, length=50), y=seq(-3*w0, 3*w0, length=50),
                             z=10))
res <- adply(xyz, 1, gaussian_field, lambda=500, maxEval = 200, 
             w0=w0, ni=1.5+0i, alpha=70*pi/180, .progress="text")
xyz <- data.frame(xyz, field=res[[2]])


p <- 
  ggplot(xyz, aes(x/1e3, y/1e3, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x="x / um", y="y / um", fill=expression("|E|"^2)) + 
  coord_fixed()

p
# 
## 3D map
# xyz <- as.matrix(expand.grid(x=seq(-20e3, 20e3, length=40), y=seq(-20e3, 20e3, length=40),
#                              z=seq(0, 15e3, length=40)))
# res <- adply(xyz, 1, field, lambda=500, w0=1e4, ni=1.5+0i, alpha=60*pi/180, .progress="text")
# 
# #save(res,file="3D.rda")
# #load(file="3D.rda")
# m <- array(res[[2]], c(40,40,40))
# 
# library(rgl)
# open3d()
# sprites3d( xyz[,1], xyz[,2], xyz[,3] - mean(xyz[,3]), radius=1000, color="blue",
#           lit=FALSE, alpha=scales::rescale(res[[2]]))
# #aspect3d(1,1,0.1)
