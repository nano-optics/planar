# ---- setup -----
library(planar)
opts_chunk$set(fig.path="tamm/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))
library(reshape2)
# prettier palette
palette(palette_tamm)

scale_colour_discrete <- function(...) 
  scale_colour_brewer(..., palette="Set1")

scale_fill_discrete <- function(...) 
  scale_fill_manual(..., values=palette())

# ---- structure -----



tamm <- tamm_stack_ir(pairs=30, dm=28, metal = "epsAu", dx1 = -20,
                      n1=3.7, n2=3.0, nleft=1, dleft=1000, position = "before")
p <- autoplot(tamm)
# ---- ff -----
ff <- simulate_ff(s=tamm, wavelength=seq(800, 1100))
head(ff)
mff <- melt(ff, meas=c("R","T","A"))
p1 <- 
  ggplot(mff, aes(wavelength, value, colour=variable))+
  geom_line() +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="")

optimum <- subset(ff, A == max(A))
optimum

# ---- nf -----
nf <- simulate_nf(s=tamm, wavelength=optimum$wavelength)
head(nf)

p2 <- 
  p + geom_line(aes(x, I), data=nf) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0)) +
  ggtitle(sprintf("wavelength: %.2f", optimum$wavelength))

library(gridExtra)

grid.arrange(p1, p2)

ggsave("scheme.pdf", arrangeGrob(p1, p2))

library(plyr)
model <- function(dx=0){
  tamm <- tamm_stack_ir(lambda0 = 1000, pairs=30, dm=30, metal = "epsAu", dx1 = dx,
                        n1=3.7, n2=3.0, nleft=1, dleft=1000, position = "before")
  ff <- simulate_ff(s=tamm, wavelength=seq(950, 1100, length=500))
 melt(ff, meas=c("A"))
}

pm <- data.frame(dx=seq(-30,-10, length=10))
all <- mdply(pm, model)

# ggplot(all, aes(1e7/885 - 1e7/wavelength, value, colour=dx, group=dx))+
p3 <- 
  ggplot(all, aes(wavelength, value, colour=dx, group=dx))+
  geom_line(lwd=1.2) +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  scale_x_continuous( expand=c(0,0)) +
  scale_colour_distiller(palette = "Spectral")+
  theme_grey() + guides(colour="none") +
  labs(x = "Wavelength /nm", y="Absorbance")
  
  # ggsave("tuning.pdf")

ggsave("scheme.pdf", width=12, height=6, arrangeGrob(p1, p2 +coord_cartesian(xlim=c(500, 3000)), p3, widths=c(1.5,2),
                                 layout_matrix = matrix(c(1,2,3,3),ncol=2,nrow=2)))
