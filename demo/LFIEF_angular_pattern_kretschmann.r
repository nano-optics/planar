## radiation pattern of a dipole near a dielectric/(metal)/dielectric interface
## parallel and perpendicular orientations, p- polarisation

library(planar)
require(reshape2)
library(ggplot2)

## from left to right
## incident water | metal | glass

params <- list(lambda=632.8,
               nPrism = 1.5,
               nWater = 1.0,
               metal = epsAg(632.8)$epsilon,
               theta = seq(-90, 90, length=300) * pi/180)

back <- list(lambda = params$lambda,
             theta = params$theta, 
             epsilon = list(params$nPrism^2, params$metal, params$nWater^2),
             thickness = c(0, 50, 0),
             d = 1,
             polarisation = 'p')


front <- invert_stack(back)

## actual simulation
M.front <- do.call(multilayer, front)
M.back <- do.call(multilayer, back)

theta <- c(params$theta + pi/2, params$theta + 3*pi/2) * 180 / pi
## look at field in first and last media
last <- length(M.back$Mr.par)
intensity.par = c(M.front$Ml.par[[1]] , M.back$Mr.par[[last]])
intensity.perp = c(M.front$Ml.perp[[1]] , M.back$Mr.perp[[last]])

combined <- data.frame(theta=theta, parallel=intensity.par,
                       perpendicular=intensity.perp,
                       side = gl(2, length(M.front$q)))

## basic plot
qplot(theta, perpendicular, colour=side, data=combined, geom="line")

## polar plot with two variables
m <- melt(combined, id=c("theta", "side"))

p <- 
ggplot(m, aes(theta, value, group=variable:side)) +
  geom_polygon(aes(fill=side), alpha=0.5)  +
  geom_path(aes(linetype=variable), alpha=0.5)  +
  scale_x_continuous(breaks = seq(0, 360, by = 45),
                     limits = c(0, 360), expand = c(0, 0)) +
  scale_y_log10(lim=c(1e-4, 200)) +
  coord_polar() +
  annotate("segment", x=360, xend=360, y=0, yend=1, colour="orange", size=2) +
  annotate("segment", x=180, xend=180, y=0, yend=1, colour="orange", size=2) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(x = "angle / degrees", y = "LFIEF", fill = "side", linetype = "dipole")

p
