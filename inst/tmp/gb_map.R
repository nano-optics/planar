library(planar)
library(plyr)
library(ggplot2)
simulation_bottom <- function(w0=1e4, angle=0,
                              wavelength=632.8,
                              epsilon=c(1.5^2, epsAg(wavelength)$epsilon, 1.0^2),
                              thickness = c(0, 50,0)){
  angle <- angle*pi/180
  
  xyz <- as.matrix(expand.grid(x=seq(-2*w0-10*wavelength*sin(angle),
                                     2*w0+50*wavelength, length=200), 
                               y=0,
                               z=seq(-10*wavelength, 2*wavelength, length=200)))
  
  res <- gaussian_near_field2(xyz, wavelength=wavelength,
                              epsilon=unlist(epsilon), thickness=thickness,
                              w0=w0, alpha=alpha, maxEval=500)
  
  data.frame(xyz, field=res)
}


simulation_top <- function(w0=1e4, angle=0,
                       wavelength=632.8,
                       epsilon=c(1.5^2, epsAg(wavelength)$epsilon, 1.0^2),
                       thickness = c(0, 50,0)){
  angle <- angle*pi/180
  
  xyz <- as.matrix(expand.grid(x=seq(-2*w0,
                                     2*w0+50*wavelength, length=200), 
                               y=0,
                               z=seq(0, wavelength, length=200)))
  
  res <- gaussian_near_field2(xyz, wavelength=wavelength,
                              epsilon=unlist(epsilon), thickness=thickness,
                              w0=w0, alpha=alpha, maxEval=500)
  
  data.frame(xyz, field=res)
}

params <- expand.grid(w0=c(5)*1e3,
                      angle=seq(42, 46, by=4))

params <- expand.grid(w0=c(2)*1e3,
                      angle=c(45))

bottom <- mdply(params, simulation_bottom, .progress="none")
top <- mdply(params, simulation_top, .progress="none")


# xyz <- simulation(w0=1e3, angle=45)

p <- ggplot(bottom, aes(x, z, fill=field))+
#   facet_grid(angle~w0, scales="free")+
  geom_raster(interpolate=TRUE) +
  geom_raster(data=top, interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /nm"), fill=expression("|E|"^2), 
       y=expression("z /nm")) +
  coord_equal()

p

# ggsave("gb_2um_45deg.png", p, width=20, height=4)
# 
# 
# plot_one <- function(d)
#   ggplot(d, aes(x/w0, z, fill=field))+
#   geom_raster() +
# #   geom_tile(aes(height=z)) +
#   scale_x_continuous(expand=c(0,0))+
#   scale_y_continuous(expand=c(0,0), lim=c(-1000, 500)) +
#   labs(x=expression("x /nm"), fill=expression("|E|"^2), 
#        y=expression("z /nm")) + guides(fill="none")
#   
# lp <- dlply(all, c("w0", "angle"), plot_one)
# library(gtable)
# lg <- llply(lp, ggplotGrob)
# lg <- llply(lg, gtable_filter, "panel")
# 
# do.call(grid.arrange, c(lg, nrow=3, as.table=TRUE))
# 
# 
