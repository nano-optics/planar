
## ----load, echo=FALSE,results='hide'-----------------------------------------
library(knitr)
library(ggplot2)
library(rgl)
knit_hooks$set(rgl = function(before, options, envir) {
  # if a device was opened before this chunk, close it
  if (before && rgl.cur() > 0) rgl.close()
  hook_rgl(before, options, envir)
})
opts_chunk$set(fig.path="gaussianbeam/", cache=TRUE, cache.path="gaussian/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----setup, results='hide', echo=FALSE----------------------------------------
library(planar)
library(ggplot2)
library(plyr)


## ----depth---------------------------------------------------------------
xyz <- as.matrix(expand.grid(x=seq(0,1e4,length=20), y=0, z=seq(1, 500, length=50)))
res <- adply(xyz, 1, gaussian_near_field, d=50, nl=sqrt(epsAu(633)$epsilon),
             wavelength=633, w0=5e3, ni=1.5, alpha=45, maxEval=1000)
xyz <- data.frame(xyz, field=res[[2]])
ggplot(xyz, aes(z, field, colour=x, group=x))+
  geom_line() + scale_y_log10()+
  labs(x=expression("z /nm"), y=expression(log("|E|"^2)), 
       colour=expression("x /nm")) 



## ----waist---------------------------------------------------------------
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
  xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0+5000,length=100), y=0, z=seq(1)))
  res <- adply(xyz, 1, gaussian_near_field, d=50, nl=sqrt(metal),
               wavelength=wavelength, w0=w0, ni=1.5, alpha=spp, maxEval=2000)
  data.frame(xyz, field=res[[2]])
}

params <- data.frame(w0=c(5, 10, 50, 1e2, 500, 1e3, 1e4 ))
all <- mdply(params, simulation)

ggplot(all, aes(x/w0/1000, field, group=w0, colour=factor(w0)))+
  geom_line() + 
  geom_vline(aes(x=0,y=NULL),lty=2) +
  geom_hline(aes(x=0,yintercept=maxi),lty=3) +
  annotate("text", label="plane-wave", y=maxi, x=-2.5, vjust=1, fontface="italic") +
  labs(x=expression(x/w[0]), y=expression("|E|"^2), 
       colour=expression(w[0]/mu*m)) +
  coord_cartesian(xlim=c(-5,5)) + theme_minimal()+
  guides(colour=guide_legend(reverse=TRUE)) +
  theme(panel.background=element_rect(fill=NA))


## ----GH------------------------------------------------------------------
w0 <- 1000
xyz <- as.matrix(expand.grid(x=seq(-5*w0, 15*w0,length=100), y=0, z=seq(0, 500, by=50)))
res <- adply(xyz, 1, gaussian_near_field, d=50, nl=sqrt(epsAu(633)$epsilon),
             wavelength=633, w0=w0, ni=1.5, alpha=44.5*pi/180, maxEval=2000)
xyz <- data.frame(xyz, field=res[[2]])

ggplot(xyz, aes(x, field, colour=z, group=z))+
  geom_line() +
  geom_vline(aes(x=0,y=NULL),lty=2) +
  labs(x=expression("x /nm"), y=expression("|E|"^2), 
       colour=expression("z /nm")) 


## ----map-----------------------------------------------------------------
w0 <- 2000
xyz <- as.matrix(expand.grid(x=seq(-3*w0, 5*w0, length=50), 
                             y=seq(-3*w0, 3*w0, length=50),
                             z=10))
res <- adply(xyz, 1, gaussian_near_field, d=50, nl=sqrt(epsAu(633)$epsilon),
             wavelength=633, w0=w0, ni=1.5, alpha=44.5*pi/180, maxEval=500)
xyz <- data.frame(xyz, field=res[[2]])

  ggplot(xyz, aes(x/1e3, y/1e3, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /nm"), fill=expression("|E|"^2), 
       y=expression("y /nm")) +
  coord_fixed()



## ----volume, rgl=TRUE----------------------------------------------------
w0 <- 2000
x <- seq(-5*w0, 15*w0, length=100)
z <- seq(2, 1000, length=100)
xyz <- as.matrix(expand.grid(x=x, 
                             y=0, 
                             z=z))
                             
res <- adply(xyz, 1, gaussian_near_field, d=50, nl=sqrt(epsAu(633)$epsilon),
             wavelength=633, w0=w0, ni=1.5, alpha=44.5*pi/180, maxEval=500)
xyz <- data.frame(xyz, field=res[[2]])

col <- heat.colors(100)[cut(xyz$field, 100)]

rgl.surface(x=x, y=matrix(xyz$field/max(xyz$field)*1e4, ncol=100), 
            z=z*10, color=col, back="lines", lit=FALSE, alpha=0.9)


