## radiation pattern of a dipole in a vacuum
## parallel and perpendicular orientations, p- polarisation

library(planar)
require(reshape2)
library(ggplot2)

## from left to right
## incident air | metal | glass

back <- list(lambda=632.8,
             theta = seq(0, 180, length=100)* pi/180, 
             epsilon = list(1.0^2, 1.0^2),
             thickness = c(0, 0),
             d = 1,
             polarisation = 'p')


front <- invert_stack(back)

## simulation
M.front <- do.call(multilayer, front)
M.back <- do.call(multilayer, back)

theta <- c(back$theta , 2*pi-rev(back$theta) ) * 180 / pi
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

mylabels <- c(-90,-45,0,45,90, 45,0,-45)
p <- 
  ggplot(m, aes(theta, value)) +
  geom_polygon(aes(fill=variable,colour=variable), alpha=0.5)  +
  scale_x_continuous(breaks = seq(0, 360-45, by = 45), labels=mylabels,
                     limits = c(0, 360), expand = c(0, 0)) +
  coord_polar(start=0) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "angle / degrees", y = "LFIEF", fill = "side",linetype="Orientation") +
  guides(fill="none",colour="none")

p 
