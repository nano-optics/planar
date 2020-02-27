
## ----load, echo=FALSE,results='hide'-----------------------------------------
library(knitr)
library(ggplot2)

library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----setup, results='hide', echo=FALSE----------------------------------------
library(planar)
library(ggplot2)
library(plyr)


## ----pw---------------------------------------------------------------
struct <- list(wavelength=632.8,
               thickness=c(0,50,0),
               epsilon=c(1.5^2, epsAu(632.8)$epsilon, 1.0^2))

w0 <- 800
xyz <- as.matrix(expand.grid(x=seq(-4*w0, 6*w0, length=300), 
                             y=0,
                             z=seq(-4*w0, w0, length=300)))

# res <- gaussian_near_field_ml(xyz, epsilon=struct$epsilon,
#                               wavelength=struct$wavelength, 
#                               thickness=struct$thickness,
#                               w0=w0, alpha=spp, maxEval=1000)
# 
# m <- data.frame(xyz, field=res)

# saveRDS(m,file='spp.rds')
m <- readRDS('spp.rds')

p <- ggplot(m, aes(x/1e3, z/1e3, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))`+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /nm"), fill=expression("|E|"^2), 
       y=expression("y /nm")) +
  # coord_fixed() +
  scale_fill_distiller(palette = 5, direction = 1)+
  guides(fill='none')+
  theme_minimal()+
  theme()

p

# ggsave('spp.png', p, width=10, height=5, dpi=300)


# library(rayshader)
# 
# plot_gg(p)

# render_movie('test.mp4')

