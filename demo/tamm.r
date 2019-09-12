# ---- setup -----
library(planar)
library(reshape2)
opts_chunk$set(fig.path="tamm/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))

# prettier palette
palette(palette_tamm)

scale_colour_discrete <- function(...) 
  scale_colour_brewer(..., palette="Set1")

scale_fill_discrete <- function(...) 
  scale_fill_manual(..., values=palette())

# ---- structure -----

tamm <- tamm_stack_ir(pairs=10, dm=50)
(pp <- autoplot(tamm))

# ---- ff -----
ff <- simulate_ff(s=tamm, wavelength=seq(600, 1200), angle = 0, polarisation = 'p')
head(ff)
mff <- melt(ff, meas=c("R","T","A"))
p0 <- ggplot(mff, aes(wavelength, value, colour=variable))+
  geom_line() +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="")

p0

optimum <- subset(ff, A == max(A))
optimum

# ---- nf -----
nf <- simulate_nf(s=tamm, wavelength=optimum$wavelength, angle = 0, polarisation = 'p')
head(nf)

pp + geom_line(aes(x, I), data=nf) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))

w0 <- 6000
xyz <- as.matrix(expand.grid(x=seq(-3*w0, 3*w0, length=50), 
                             y=0,
                             z=seq(-3000, 2000, length=200)))


# xyz <- as.matrix(expand.grid(x=seq(-3*w0, 3*w0, length=30), 
#                              y=0,
#                              z=seq(-3000, 3000, length=30)))
# 
# # 
# xyz <- as.matrix(expand.grid(x=0,
#                              y=0,
#                              z=seq(-3000, 3000, length=200)))

leps <- epsilon_dispersion(tamm[["epsilon"]], optimum$wavelength)
res <- gaussian_near_field_ml(xyz, epsilon=unlist(leps),
                              wavelength=optimum$wavelength,
                              thickness=tamm$thickness, psi = 0,tol = 1e-03,
                              w0=w0, alpha=0, maxEval=200, progress = TRUE)

m <- data.frame(xyz, field=res)

# ggplot(m, aes(z, field)) + geom_line()


# saveRDS(m,file='tamm.rds')


p <- ggplot(m, aes(x/1e3, z/1e3, fill=field))+
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression("x /µm"), fill=expression("|E|"^2), 
       y=expression("y /µm")) +
  # coord_fixed() +
  scale_fill_distiller(palette = 5, direction = 1)+
  guides(fill='none') +
  theme_minimal() +
  theme()

p
library(rayshader)

plot_gg(p)

