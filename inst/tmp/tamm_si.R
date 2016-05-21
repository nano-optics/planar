# ---- setup -----
library(planar)
opts_chunk$set(warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)
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

tamm_stack_si <- function (lambda0 = 450, n1 = 1.8, n2 = 2.8, d1 = lambda0/4/n1, 
          d2 = lambda0/4/n2, N = 2 * pairs, pairs = 4, dx1 = 0, dx2 = 0, 
          dm = 50, metal = "epsSi", position = "after", incidence = "right", 
          nleft = n2, nright = 1, ...) 
{
  tamm_stack(lambda0 = lambda0, n1 = n1, n2 = n2, d1 = d1, 
             d2 = d2, N = N, pairs = pairs, dx1 = dx1, dx2 = dx2, 
             dm = dm, metal = metal, position = position, incidence = incidence, 
             nleft = nleft, nright = nright, ...)
}

tamm <- tamm_stack_si(pairs=11, metal = "epsAg", n1 = 2, n2=3, dm=500, dx1=0, incidence="left", position="after")
(p <- autoplot(tamm))

# ---- ff -----
ff <- simulate_ff(s=tamm, wavelength=seq(300, 600, length=1e3))
head(ff)
mff <- melt(ff, meas=c("R","T","A"))
ggplot(mff, aes(wavelength, value, colour=variable))+
  geom_line() +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="")



search <- function(d, metal = "epsAg",n1 = 2, n2=3, pairs=10, wavelength=seq(350, 450, length=1e3), ...){
  
  tamm <- tamm_stack_si(pairs=pairs, metal = metal, 
                        n1 = n1, n2=n2, dm=d, incidence="left", position="before", ...)

  melt(simulate_ff(s=tamm, wavelength=wavelength), id="wavelength", meas=c("R", "A"))
}

pa <- expand.grid(d=seq(30, 100, length=10))
all <- mdply(pa, search, wavelength=seq(380, 420, length=1e3))

ggplot(all, aes(wavelength, value, colour=d, group=d)) + 
  facet_wrap(~variable) +
  geom_line()


pa <- expand.grid(d=seq(0, 10, length=10))
all <- mdply(pa, search, n1=3.7, n2=3, pairs=30, wavelength=seq(350, 500, length=1e3), metal="epsSi")

ggplot(all, aes(wavelength, value, colour=d, group=d)) + 
  facet_wrap(~variable) +
  geom_line() +
  scale_y_continuous(expand=c(0,0), lim=c(0,1))



optimum <- subset(ff, A == max(A))
optimum

# ---- nf -----
nf <- simulate_nf(s=tamm, wavelength=optimum$wavelength)
head(nf)

p + geom_line(aes(x, I), data=nf) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))

