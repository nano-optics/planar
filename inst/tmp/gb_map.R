library(planar)
library(plyr)
library(ggplot2)

simulation_bottom <- function(w0=1e4, angle=0,
                              wavelength=632.8,
                              epsilon=c(1.5^2, epsAg(wavelength)$epsilon, 1.0^2),
                              thickness = c(0, 50,0)){
  angle <- angle*pi/180
  
  xyz <- as.matrix(expand.grid(x=seq(-2*w0-10*wavelength*sin(angle),
                                     2*w0+50*wavelength, length=100), 
                               y=0,
                               z=seq(-15*wavelength, 2*wavelength, length=100)))
  
  res <- gaussian_near_field_ml(xyz, wavelength=wavelength,
                              epsilon=unlist(epsilon), thickness=thickness,
                              w0=w0, alpha=angle, maxEval=500)
  
  data.frame(xyz, field=res)
}


simulation_top <- function(w0=1e4, angle=0,
                       wavelength=632.8,
                       epsilon=c(1.5^2, epsAg(wavelength)$epsilon, 1.0^2),
                       thickness = c(0, 50,0)){
  angle <- angle*pi/180
  
  xyz <- as.matrix(expand.grid(x=seq(-2*w0,
                                     2*w0+50*wavelength, length=100), 
                               y=0,
                               z=seq(0, wavelength, length=100)))
  
  res <- gaussian_near_field_ml(xyz, wavelength=wavelength,
                              epsilon=unlist(epsilon), thickness=thickness,
                              w0=w0, alpha=angle, maxEval=500)

  
  
  data.frame(xyz, field=res)
}

params <- expand.grid(w0=c(5)*1e3,
                      angle=42)


bottom <- mdply(params, simulation_bottom, .progress="text")
top <- mdply(params, simulation_top, .progress="text")

both <- rbind(bottom, top)


p <- ggplot(top, aes(x, z, fill=field))+
  # facet_grid(angle~w0, scales="free")+
  geom_raster(data = top, interpolate=TRUE) +
  geom_raster(data = bottom, interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /nm"), fill=expression("|E|"^2),
       y=expression("z /nm")) +
  coord_equal()

p