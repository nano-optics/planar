## radiation pattern of a dipole near a dielectric/(metal)/dielectric interface
## parallel and perpendicular orientations, p- polarisation

library(planar)
require(reshape2)
library(ggplot2)

## from left to right
## incident air | metal | glass

back <- list(wavelength=632.8,
             angle = seq(-90, 90, length=1e4)[-c(1, 1e4)] * pi/180, 
             epsilon = list(1.52^2, epsAg(632.8)$epsilon, 1.0^2),
             thickness = c(0, 50, 0),
             d = 1,
             polarisation = 'p')


front <- invert_stack(back)

## simulation
M.front <- do.call(multilayer, front)
M.back <- do.call(multilayer, back)

theta <- c(back$angle + pi/2, back$angle + 3*pi/2) * 180 / pi
## look at field in first and last media
last <- length(M.back$Mr.par)
intensity.par = c(M.front$Ml.par[[1]] , M.back$Mr.par[[last]])
intensity.perp = c(M.front$Ml.perp[[1]] , M.back$Mr.perp[[last]])

combined <- data.frame(angle=theta, parallel=intensity.par,
                       perpendicular=intensity.perp,
                       side = gl(2, length(M.front$q)))

## basic plot
qplot(angle, perpendicular, colour=side, data=combined, geom="line")

## polar plot with two variables
m <- melt(combined, id=c("angle", "side"))

mylabels <- c(-90,-45,0,45,90, 45,0,-45)
p <- 
  ggplot(m, aes(angle, value, group=interaction(side,variable))) +
  annotate("rect", xmin=180,xmax=360,ymin=1e-2,ymax=100,fill="grey95",alpha=0.5) +
  annotate("rect", xmin=0,xmax=360,ymin=1e-2,ymax=1,fill=NA,colour="grey50",lty="dashed") +
  geom_polygon(aes(fill=side), alpha=0.5,colour=NA)  +
  geom_line(aes(linetype=variable,colour=side), alpha=0.5)  +
  scale_x_continuous(breaks = seq(0, 360-45, by = 45), labels=mylabels,
                     limits = c(0, 360), expand = c(0, 0)) +
  scale_y_log10(lim=c(1e-2, 200)) +
  coord_polar(start=3*pi/2) +
  annotate("segment", x=360, xend=360, y=0, yend=2, colour="orange", size=2) +
  annotate("segment", x=180, xend=180, y=0, yend=2, colour="orange", size=2) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "angle / degrees", y = "LFIEF", fill = "side",linetype="Orientation") +
  guides(fill="none",colour="none")

p
