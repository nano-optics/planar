# ---- setup -----
library(planar)
library(ggplot2)
library(reshape2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))

# prettier palette
palette(palette_tamm_porous)

scale_colour_discrete <- function(...) 
  scale_colour_brewer(..., palette="Set1")

scale_fill_discrete <- function(...) 
  scale_fill_manual(..., values=palette())

wavenumber <- function(nm){
  1e7 / nm
}

# ---- structure -----


dbr <- tamm_stack_porous(N=8, dm=0, incidence = 'left', nleft = 1, nright=1.5)
tamm <- tamm_stack_porous(N=8, dm=12, incidence = 'left', nleft = 1, nright=1.5)
(pp <- autoplot(tamm))

# ---- ff -----
ff_ref <- simulate_ff(s=dbr, wavelength=seq(400, 900,length=1000), angle = 0, polarisation = 'p')
ff <- simulate_ff(s=tamm, wavelength=seq(400, 900,length=1000), angle = 0, polarisation = 'p')
head(ff)
mff <- melt(ff, meas=c("R","T","A"))
mff2 <- melt(ff_ref, meas=c("R","T","A"))
p0 <- ggplot(mff, aes(wavelength, value, colour=variable))+
  geom_line() +
  geom_line(data=mff2,lty=2) +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="")

ggsave('ff.pdf',p0,width=5, height=3)

optimum <- subset(ff, R == min(R))
optimum

# ---- nf -----
nf <- simulate_nf(s=tamm, wavelength=optimum$wavelength, 
                  angle = 0, polarisation = 'p', dmax = 100)
head(nf)

pp2 <- pp + geom_line(aes(x, I), data=nf) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))


ggsave('nf.pdf',pp2,width=5, height=3)


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

