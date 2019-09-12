# ---- setup -----
library(planar)
library(ggplot2)
library(reshape2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))

# prettier palette
palette(palette_tamm)

scale_colour_discrete <- function(...) 
  scale_colour_brewer(..., palette="Set1")

scale_fill_discrete <- function(...) 
  scale_fill_manual(..., values=palette())

wavenumber <- function(nm){
  1e7 / nm
}

# ---- structure -----
sic <- function(wavelength, eps_inf = 6.7, nu_lo = 969, nu_to = 793, gamma = 4.76){
 nu <- wavenumber(wavelength)
 data.frame(wavelength = wavelength,
            epsilon = eps_inf * (nu_lo^2 - nu^2 - 1i*gamma*nu) / (nu_to^2 - nu^2 - 1i*gamma*nu))
}

nu0 <- 870
1e7 / 870
tamm <- tamm_stack_ir(pairs = 6, dm=1e6, metal = "sic",nleft = 1, nright = 1, 
                      lambda0 = wavenumber(870),
                      n2 = 2.38, n1 = 1.59, d2 = 1270, d1 = 1860, dx2 = -600)
(p <- autoplot(tamm))

wavelength=seq(7e3, 6e4)
test <- sic(wavelength)[['epsilon']]
matplot(wavenumber(wavelength), cbind(Re(test), Im(test)), t='l',lty=1,ylim=c(-5,10))

# ---- ff -----
ff <- simulate_ff(s=tamm, wavelength = wavenumber(seq(500, 1500)))

substrate <- planar::layer_stack(epsilon = "sic", thickness = 1e6)
ff <- simulate_ff(s=substrate, wavelength = wavenumber(seq(500, 1500)))
# ff <- simulate_ff(s=tamm, wavelength = seq(500, 1500))



head(ff)
mff <- melt(ff, meas=c("R","T","A"))
ggplot(mff, aes(wavenumber(wavelength), value, colour=variable))+
  geom_line() +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  scale_x_reverse() +
  labs(x = expression("wavenumber /"*cm^-1), y="")

ggplot(mff, aes(wavelength, value, colour=variable))+
  geom_line() +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  # scale_x_reverse() +
  labs(x = expression("wavelength /"*nm), y="")



optimum <- subset(ff, A == max(A))
optimum

# ---- nf -----
nf <- simulate_nf(s=tamm, wavelength=optimum$wavelength)
head(nf)

p + geom_line(aes(x, I), data=nf) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))

