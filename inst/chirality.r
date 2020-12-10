## ----load, echo=FALSE,results='hide'-----------------------------------------
library(knitr)
library(ggplot2)
opts_chunk$set(fig.path="internalfieldcomparison/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----setup, results='hide'----------------------------------------------------
library(planar)
library(ggplot2)
require(reshape2)
require(plyr)


## ----simulation----------------------------------------------------------

s <- tamm_stack(lambda0=705, pairs=9, incidence = "left",
                nH = 1.7, nL = 1.3, dm = 200,
                nleft = 1.5, nright = 1.0)

s2 <- embed_stack(layer_stack(epsilon = 1.0^2),dleft = 1e3,
                nleft = 1.0, nright = 4)

plot(s2)
# plot(s)
ff <- simulate_ff(s=s)
ggplot(ff, aes(wavelength, R)) + geom_line() + scale_y_continuous(lim=c(0,1))

peak <- ff$wavelength[which.min(ff$R)]

res <- internal_field(wavelength=peak, angle=0, psi=0,
                       thickness = s$thickness,
                       dmax=100,  res=1e3,
                       epsilon=unlist(epsilon_dispersion(s$epsilon, wavelength = peak)),
                       field = TRUE)

limits2 <- ddply(res, .(id), summarize,  material=unique(material),
                 xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-Inf, 
                 ymax=Inf)
library(dplyr)
library(tidyr)

mm <- res %>% 
  mutate(E2 = Mod(Ex)^2 + Mod(Ey)^2 + Mod(Ez)^2,
         H2 = Mod(Hx)^2 + Mod(Hy)^2 + Mod(Hz)^2,
         C = Im(Conj(Ex)*Hx + Conj(Ey)*Hy + Conj(Ez) * Hz),
         g = C/E2) %>% 
  select(x,H2,E2,C,g) %>% 
  pivot_longer(c('H2','E2','C','g'))
  
  # cols <- RColorBrewer::brewer.pal(5,'Set1')

ggplot(s) + 
  geom_rect(data=limits2, aes(xmin=xmin, xmax=xmax, 
                              ymin=-Inf, ymax=Inf, 
                              fill=material),alpha=0.4)+
  # geom_line(aes(x, sum), lwd=1.2) +
  geom_abline(intercept=0,slope=0,lty=3)+
  geom_line(data=mm, aes(x, y=value, 
                           colour=name,group=name,linetype=name),alpha=1) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank()) +
  guides(fill = guide_legend(override.aes=list(size=1)),
         alpha='none')+
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous(expression("|E|"^2~",|H|"^2~",C")) +
  scale_fill_brewer(palette="Pastel1")+
  scale_colour_brewer(palette="Set1")+
  scale_linetype_manual(values = c('E2'=2,'H2'=2,'C'=1,'g'=1))+
  scale_alpha(range = c(0.6,0.6)) +
  labs(colour='variable',linetype='variable')


## ---

res2 <- internal_field(wavelength=peak, angle=0, psi=0,
                      thickness = s2$thickness,
                      dmax=300,  res=1e3,
                      epsilon=unlist(epsilon_dispersion(s2$epsilon, wavelength = peak)),
                      field = TRUE)

limits2 <- ddply(res2, .(id), summarize,  material=unique(material),
                 xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-Inf, 
                 ymax=Inf)
library(dplyr)
library(tidyr)

mm <- res2 %>% 
  mutate(E2 = Mod(Ex)^2 + Mod(Ey)^2 + Mod(Ez)^2,
         H2 = Mod(Hx)^2 + Mod(Hy)^2 + Mod(Hz)^2,
         C = Im(Conj(Ex)*Hx + Conj(Ey)*Hy + Conj(Ez) * Hz),
         g = C/E2) %>% 
  select(x,H2,E2,C,g) %>% 
  pivot_longer(c('H2','E2','C','g'))

# cols <- RColorBrewer::brewer.pal(5,'Set1')

ggplot(s2) + 
  geom_rect(data=limits2, aes(xmin=xmin, xmax=xmax, 
                              ymin=-Inf, ymax=Inf, 
                              fill=material),alpha=0.4)+
  # geom_line(aes(x, sum), lwd=1.2) +
  geom_abline(intercept=0,slope=0,lty=3)+
  geom_line(data=mm, aes(x, y=value, 
                         colour=name,group=name,linetype=name),alpha=1) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  guides(fill = guide_legend(override.aes=list(size=1)),
         alpha='none')+
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous(expression("|E|"^2~",|H|"^2~",C")) +
  scale_fill_brewer(palette="Pastel1")+
  scale_colour_brewer(palette="Set1")+
  scale_linetype_manual(values = c('E2'=2,'H2'=2,'C'=1,'g'=1))+
  scale_alpha(range = c(0.6,0.6)) +
  labs(colour='variable',linetype='variable')






