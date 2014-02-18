library(planar)
library(ggplot2)
library(plyr)

theme_set(theme_minimal() + theme(panel.background=element_rect(fill=NA)))

## line example
xyz <- as.matrix(expand.grid(x=seq(0,5e2,length=3), 
                             y=10, z=seq(-1000, 1000, length=100)))

test <- collection_ml(xyz, omega=c(0,90)*pi/180, wavelength=600,
                      maxEval=500.0)


xyz <- data.frame(xyz, I=drop(test))
ggplot(xyz, aes(z, I, colour=x, group=x))+
  geom_line()


xyz <- as.matrix(expand.grid(x=2*seq(-200,200,length=30), 
                             y=2*seq(-200,200,length=30), z=1))

test <- collection_ml(xyz, omega=c(0,50)*pi/180, maxEval=500.0)


xyz <- data.frame(xyz, I=drop(test))
ggplot(xyz, aes(x, y, fill= I))+
  geom_raster()


xyz <- as.matrix(expand.grid(x=2*seq(-200,200,length=30), 
                             z=2*seq(-500,300,length=30), y=0))

test <- collection_ml(xyz, omega=c(0,80)*pi/180, 
                      maxEval=500.0, psi=0)


xyz <- data.frame(xyz, I=drop(test))
ggplot(xyz, aes(x, z, fill= I))+
  geom_raster()


