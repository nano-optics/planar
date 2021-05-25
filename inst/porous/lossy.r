# ---- setup -----
library(planar)
library(ggplot2)
library(reshape2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))

# prettier palette
# dput(palette_tamm_porous)
mypal <- c("yellowgreen","goldenrod1", "wheat", "grey98", "#B2DF8A", "#D8E8F0", "#79C171")
# palette(mypal)

scale_colour_discrete <- function(...) 
  scale_colour_brewer(..., palette="Set1")

scale_fill_discrete <- function(...) {
  scale_fill_brewer(..., palette="Pastel1")
  }
  # scale_fill_manual(..., values=palette())

wavenumber <- function(nm){
  1e7 / nm
}

rawTi <- read.csv('Ti_Johnson.csv')
rawTi$epsilon <- complex(real=rawTi$n, imaginary = rawTi$k)^2
rawCr <- read.csv('Cr_Johnson.csv')
rawCr$epsilon <- complex(real=rawCr$n, imaginary = rawCr$k)^2
str(rawTi)
str(rawCr)

epsTi <- function(wavelength){
  
  re <- spline(x = rawTi$wl*1e3, y = Re(rawTi$epsilon), xout = wavelength)
  im <- spline(x = rawTi$wl*1e3, y = Im(rawTi$epsilon), xout = wavelength)
  
  data.frame(wavelength = wavelength, epsilon = re$y + 1i*im$y)
  
}


epsCr <- function(wavelength){
  
  re <- spline(x = rawCr$wl*1e3, y = Re(rawCr$epsilon), xout = wavelength)
  im <- spline(x = rawCr$wl*1e3, y = Im(rawCr$epsilon), xout = wavelength)
  
  data.frame(wavelength = wavelength, epsilon = re$y + 1i*im$y)
  
}


print(epsCr(633))
print(epsTi(633))

# plot(epsCr(seq(400,800)),t='l')
# points(rawCr$wl*1e3,rawCr$epsilon)
# 
# 
# plot(epsTi(seq(400,800)),t='l')
# points(rawTi$wl*1e3,rawTi$epsilon)


# ---- structure -----

tamm_stack_seeded <- function(lambda0 = 600, 
                              n1 = 1.72, n2 = 1.28, 
                              d1 = lambda0/4/n1, d2 = lambda0/4/n2, 
                              N = 8, 
                              seed = "epsTi", dseed = 0, 
                              metal = "epsAu", dmetal = 12, 
                              ncoat = 1.75, dcoat = 20,
                              nleft = 1, nright = 1.5, 
                              dleft = 200, dright = 200, 
                              incidence = c("left", "right"), ...) {
  
  incidence <- match.arg(incidence)
  
  # dbr layers
  dbr_e <- as.list(rep(c(n1^2, n2^2), length.out = N))
  dbr_t <- rep(c(d1, d2), length.out = N)
  s_dbr <- structure(list(epsilon = dbr_e, thickness = dbr_t),
                   class = "stack")
  
  # coatings
  s_coat <- layer_stack(epsilon = ncoat^2, thickness = dcoat)
  s_metal <- layer_stack(epsilon = metal, thickness = dmetal)
  s_seed <- layer_stack(epsilon = seed, thickness = dseed)
  
  # full stack
  s <- c(s_coat, s_metal, s_seed, s_dbr)
  
  s <- embed_stack(s, nleft = nleft, nright = nright, 
                   dleft = dleft, dright = dright)
  
  switch(incidence, left = s, right = rev(s))
  
}

dbr <- tamm_stack_seeded(dmetal = 0, dseed = 0, dcoat = 0, incidence = 'left')
# standard tamm
tamm1 <- tamm_stack_seeded(dmetal = 12, dseed = 0, dcoat = 0, incidence = 'left')
# protective coat
tamm2 <- tamm_stack_seeded(dmetal = 12, dseed = 0, dcoat = 20, incidence = 'left')
# nefarious seeds
tamm3 <- tamm_stack_seeded(dmetal = 12, dseed = 5, seed = "epsTi",
                           dcoat = 0, incidence = 'left')
tamm4 <- tamm_stack_seeded(dmetal = 12, dseed = 5, seed = "epsCr",
                           dcoat = 0, incidence = 'left')
(pp <- autoplot(tamm4))

# grid::grid.raster(mypal)



# ---- ff -----
ff_ref <- simulate_ff(s=dbr, wavelength=seq(400, 900,length=1000), 
                      angle = 0, polarisation = 'p')

ff1 <- simulate_ff(s=tamm1, wavelength=seq(400, 900,length=1000), 
                  angle = 0, polarisation = 'p')
ff2 <- simulate_ff(s=tamm2, wavelength=seq(400, 900,length=1000), 
                  angle = 0, polarisation = 'p')
ff3 <- simulate_ff(s=tamm3, wavelength=seq(400, 900,length=1000), 
                   angle = 0, polarisation = 'p')
ff4 <- simulate_ff(s=tamm4, wavelength=seq(400, 900,length=1000), 
                   angle = 0, polarisation = 'p')

ff_all <- rbind(cbind(ff1, coat = 'none', seed = "none"),
                cbind(ff2, coat = '20',   seed = "none"),
                cbind(ff3, coat = 'none', seed = "epsTi"),
                cbind(ff4, coat = 'none', seed = "epsCr"))
head(ff)
m_ref <- melt(ff_ref, meas=c("R","T","A"))
m_all <- melt(ff_all, meas=c("R","T","A"))

p1 <- ggplot(subset(m_all, variable != 'A'), aes(wavelength, value))+
  geom_area(data=subset(m_ref, variable != 'A'),fill='grey92') +
  geom_line(aes(linetype=coat, colour=seed)) +
  facet_wrap(~variable)+
  scale_x_continuous(lim=c(400,800), expand=c(0,0)) +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="") +
  scale_linetype_manual(values=c(2,1))+
  scale_colour_manual(values=c(none='black',
                               epsCr = "#E41A1C", epsTi = "#377EB8"))+
  egg::theme_article() +
  theme(panel.spacing.x = unit(0.5,'cm'))

p1


# ggsave('seeding_o.pdf',p1,width=8, height=4)


RColorBrewer::brewer.pal(3,'Set1')


# ggsave('ff.pdf',p0,width=5, height=3)

opt1 <- subset(ff1, R == min(R))
opt2 <- subset(ff2, R == min(R))
opt3 <- subset(ff3, R == min(R))
opt1
opt2
opt3

# ---- nf -----
nf1 <- simulate_nf(s=tamm1, wavelength=694.2943, 
                  angle = 0, polarisation = 'p', dmax = 100)
nf2 <- simulate_nf(s=tamm2, wavelength=703.8038, 
                   angle = 0, polarisation = 'p', dmax = 100)
nf3 <- simulate_nf(s=tamm3, wavelength=691.2913, 
                   angle = 0, polarisation = 'p', dmax = 100)

pp1 <- autoplot(tamm1) + geom_line(aes(x, I), data=nf1) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))
pp2 <- autoplot(tamm2) + geom_line(aes(x, I), data=nf2) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))
pp3 <- autoplot(tamm3) + geom_line(aes(x, I), data=nf3) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))

library(patchwork)
pp1 + pp2 + pp3 + plot_layout(ncol=1)

# ggsave('nf_comp.pdf',pp2,width=5, height=3)


# wtf <-  lfief(wavelength = optimum$wavelength, 
#                 angle = 0, 
#                 polarisation = "p", 
#                 thickness = tamm$thickness, 
#                 dmax = 100, res = 100,  res2=1e3,
#                 epsilon = epsilon_dispersion(tamm$epsilon, wavelength = optimum$wavelength), 
#               displacement = FALSE)
# 
# 
# wtf2 <- internal_field(wavelength = optimum$wavelength, angle = 0, 
#                       thickness = tamm$thickness, 
#                       epsilon = unlist(epsilon_dispersion(tamm$epsilon, wavelength = optimum$wavelength)),
#                       psi = 0)
# 
# ggplot(wtf, aes(x, value)) +
#   geom_line(aes(group=variable))+
#   geom_line(data=wtf2, aes(x, I),col='red',lty=2)


# ---- variable_coverlayer -----

tamm_stack_covered <- function (lambda0 = 600, n1 = 1.72, n2 = 1.28, d1 = lambda0/4/n1, 
                                d2 = lambda0/4/n2, N = 8, pairs = N/2, dx1 = 0, dx2 = 0, 
                                dm = 12, metal = "epsAu", cover_thickness = 20, cover_index = 1.75,
                                position = c("before", 'after'), incidence = c('left',"right"), 
                                nleft = 1.0, nright = 1.5, dleft = 200, dright = 200, ...) {
  position <- match.arg(position)
  incidence <- match.arg(incidence)
  dbr <- dbr_stack(lambda0 = lambda0, n1 = n1, n2 = n2, d1 = d1, 
                   d2 = d2, N = N - 2, pairs = pairs - 1)
  variable <- dbr_stack(lambda0 = lambda0, n1 = n1, n2 = n2, 
                        d1 = d1 + dx1, d2 = d2 + dx2, N = 2)
  met <- layer_stack(epsilon = metal, thickness = dm)
  cover <- layer_stack(epsilon = cover_index^2, thickness = cover_thickness)
  if (N%%2 && position == "after") 
    variable <- rev(variable)
  struct <- switch(position, before = c(cover, met, variable, dbr), 
                   after = c(dbr, variable, met, cover))
  s <- embed_stack(struct, nleft = nleft, nright = nright, 
                   dleft = dleft, dright = dright)
  switch(incidence, left = s, right = rev(s))
}


palette_covered <- append(palette_tamm_porous, blues9[5], after = 5)
palette(palette_covered)

s <- tamm_stack_covered(N=8, dm=12, incidence = 'left', nleft = 1, nright=1.5)
(pp <- autoplot(s))

model <- function(cover_thickness = 20, cover_index = 1.75, ..., 
                  wavelength=seq(400, 900,length=300),
                  angle = 0, polarisation = 'p'){
  
  s <- tamm_stack_covered(cover_thickness = cover_thickness, cover_index = cover_index, ...)
  
  ff <- simulate_ff(s=s, wavelength=wavelength, angle = angle, polarisation = polarisation)
  
  melt(ff, meas=c("R","T","A"))
  
}

library(dplyr)
params <- expand.grid(cover_thickness=seq(20,100,by=20), cover_index=1.75)
params$case <- as.character(seq(1,nrow(params)))
library(tidyr)
library(purrr)

all <- params %>% pmap_df(., model, dm=12, incidence = 'left', nleft = 1, nright=1.5, .id = 'case')
results <- left_join(params, all, "case")

p <- ggplot(results %>% filter(variable !='A'), aes(wavelength, value)) +
  facet_wrap(~variable) +
  geom_area(data = mff2 %>% filter(variable !='A'), fill='grey70',alpha=0.5)+
  geom_line(aes(colour=factor(cover_thickness)), lwd=1.0) +
  geom_line(data = mff %>% filter(variable !='A'), col='black',lty=3, lwd=1.0)+
  # geom_line(data = mff2, col='black',lty=2)+
  scale_colour_hue() +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  scale_x_continuous(lim=c(400,800), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="",colour='coating /nm') +
  theme(panel.spacing.x = unit(1,'cm'))

ggsave('sweep.pdf',p,width=8, height=4)

