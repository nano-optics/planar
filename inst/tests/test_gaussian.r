library(planar)


planar:::gaussian$integrand_gb(c(0.1, 0.3), c(2.1, 3, 4), 2*pi/500, 0,pi/2 - 15*pi/180,  10, 3.2, 1.04)


f <- 1
library(cubature)
field <- function(x=1, y=1, z=1, lambda=500, alpha=pi/2 - 15*pi/180, psi=0, gamma=10, np=1.5, ns=1){

  k0 <- 2*pi/lambda
  kp <- k0*np
  
  res <- adaptIntegrate(planar:::gaussian$integrand_gb,
                        lowerLimit=c(-0.4, -0.4)*f, # some function of alpha
                        upperLimit=c(0.4, 0.4)*f, # some function of alpha
                        fDim = 6,
                        maxEval = 200,
                        r2 = c(x, y, z), kp=kp, psi=0, alpha=alpha,
                        gamma=gamma, np=np, ns=ns)$integral

  E <- complex(real = res[1:3], imag=res[4:6])
  Re(crossprod(E, Conj(E)))

}

## xyz <- as.matrix(expand.grid(x=seq(1, 2000, length=50), y=seq(1, 2000, length=50), z=4))
## res <- adply(xyz, 1, field,gamma=2,  .progress="text")
## xyz <- data.frame(xyz, field=res[[2]])
## ggplot(xyz, aes(x, y, fill=field))+
##   geom_raster(interpolate=TRUE)



xyz <- as.matrix(expand.grid(x=seq(1, 1000, length=10), y=10, z=seq(1, 1000, length=100)))
res <- adply(xyz, 1, field, lambda=500, gamma=2, np=1.5, alpha=5*pi/12,  .progress="none")
xyz <- data.frame(xyz, field=res[[2]])
ggplot(xyz, aes(z, field, colour=x, group=x))+
  geom_line()


xyz <- as.matrix(expand.grid(x=seq(-1000, 1000, length=200), y=seq(-1000, 1000, length=200),
                             z=10))
res <- adply(xyz, 1, field, lambda=500, gamma=2, np=1.5, alpha=5*pi/12, .progress="text")
xyz <- data.frame(xyz, field=res[[2]])

## xyz1 <- xyz

p <- 
ggplot(xyz, aes(x, y, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x="x / nm", y="y / nm", fill=expression("|E|"^2)) +
  theme(aspect.ratio=1)

p

## ggsave(p, file="GHxy.pdf", width=8, height=8)


## exponential <- function(gamma, s1x, s1y=s1x){

##   s1z <- sqrt(1 - s1x^2 - s1y^2+0i)
##   exp(gamma^2*(1 - s1z^(-2)))
## }


## s1 <- expand.grid(s1x=seq(-0.3, 0.3, length=100), s1y=seq(-0.3, 0.3, length=100))

## s1$f <- exponential(10, s1[, 1], s1[, 2])
