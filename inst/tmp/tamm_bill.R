setwd("~/Documents/github/planar/inst/tmp")

# ---- setup -----
library(planar)
library(reshape2)
library(plyr)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))

# prettier palette
palette(palette_tamm)

scale_colour_discrete <- function(...) 
  scale_colour_brewer(..., palette="Set1")

scale_fill_discrete <- function(...) 
  scale_fill_manual(..., values=palette())

# ---- structure -----

library(dielectric)
data("aSi")
epsSi <- function(wavelength) {
  eps <- aSi$predict(new.wavelength = wavelength)
  # eps$epsilon <- eps$epsilon + 2i
  eps
}

# dye <- readRDS("TDBC.rds")
dye <- readRDS("bill.rds")


dye_fun <- function(wavelength){
  eps_r <- predict(smooth.spline(dye$wavelength, Re(dye$epsilon)), wavelength)$y
  eps_i <- predict(smooth.spline(dye$wavelength, Im(dye$epsilon)), wavelength)$y
  data.frame(wavelength = wavelength, epsilon = eps_r + 1i*eps_i, real=eps_r, imag=eps_i)
}

test <- dye_fun(seq(300, 800))


m <- reshape2::melt(test, id="wavelength", meas=c("real", "imag"))
ggplot(m, aes(wavelength, value))+ facet_wrap(~variable, scales="free", ncol=1) + geom_line()

tamm <- tamm_stack_si(pairs=20, metal = "dye_fun", n1 = 3, n2=2, dm=55, dx1=0, incidence="left", position="before")


tamm_stack_si <- function (lambda0 = 490, 
                           n1 = 2.02, n2 = 1.46, 
                           d1 = lambda0/4/n1, d2 = lambda0/4/n2, 
                           N = 2 * pairs, pairs = 4, 
                           ds = 0, ns = n1, 
                           dm = 50, 
                           metal = "dye_fun", 
                           nleft = n2, nright = 1, 
                           dleft=500, dright=500, ...) 
{
  
  dbr <- dbr_stack(lambda0=lambda0, 
                   n1=n1, n2=n2, 
                   d1=d1, d2=d2,
                   N=N, pairs=pairs)
  
  spacer <- layer_stack(epsilon = ns^2, thickness = ds)
  met <- layer_stack(epsilon=metal, thickness=dm)
  
  struct <- c(dbr, spacer, met)
  
  s <- embed_stack(struct, nleft=nleft, nright=nright, 
                   dleft= dleft, dright=dright) 
  s
}

# 490/4/2

d <- 61
tamm <- tamm_stack_si(pairs=4, lambda0 = 510, metal = "dye_fun", dm=500, ds=d, ns=2.02)
ref <- tamm_stack_si(pairs=4, lambda0 = 510, metal = "dye_fun", dm=0, ds=d, ns=2.02)

p <- autoplot(tamm)
p 

# ---- ff -----
ff <- simulate_ff(s=tamm, wavelength=seq(400, 700, length=1e3))
fref <- simulate_ff(s=ref, wavelength=seq(400, 700, length=1e3))
head(ff)
mff <- melt(ff, meas=c("R","T","A"))
mffref <- melt(fref, meas=c("R","T","A"))

ggplot(mff, aes(wavelength, value, colour=variable))+
  geom_line() +
  geom_line(data=mffref, lty=2) +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="")


optimum <- subset(subset(ff, wavelength>450 & wavelength < 590), A == max(A))
optimum

# ---- nf -----
nf <- simulate_nf(s=tamm, wavelength=optimum$wavelength)
head(nf)

p + geom_line(aes(x, I), data=nf) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))



