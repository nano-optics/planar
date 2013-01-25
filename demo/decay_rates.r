## Reproducing Fig. 6.1, p. 304 from P&E's Principles of Surface-Enhanced Raman Spectroscopy
## a dipole is placed near a semi-infinite air/metal interface
## with orientation either parallel or perpendicular to the interface
## the total decay rates peak at the wavelength of excitation of planar SPPs epsilon=-1 at the interface (loss channel)
## the radiative decay rate in the upper medium has a trough at the wavelength where epsilon=0 (Dn=0, by continuity En=0)

library(planar)
library(ggplot2)
library(lattice)
library(gridExtra)
require(reshape2)
require(plyr)

wvl <- seq(200, 1000,by=2)*1e-3
silver <- epsAg(wvl*1e3)
gold <- epsAu(wvl*1e3)


distance <- function(d, material="silver", ...){

  material <- get(material)
  
  params <- list(d=d,
                 lambda = material$wavelength,
                 epsilon = list(incident=1.0^2, material$epsilon),
                 thickness = c(0, 0))
  
  dl <- do.call(dipole, c(params, list(Nquadrature1 = 1e3, Nquadrature2 = 5e3,
                                        Nquadrature3 = 5e3, qcut = NULL, rel.err=1e-3,
                                       show.messages = FALSE)))
  
  m <- melt(dl, id = "wavelength")
  ## print(paste(c("evaluations", dl$evaluations), collapse=" "))
  ## print(paste(c("Max rel. error", round(100*max(unlist(dl$errors)),1), "%"), collapse=" "))
  
  m$orientation <- m$variable
  
  levels(m$orientation) <- list(perpendicular="Mtot.perp",
                                perpendicular="Mrad.perp",
                                parallel="Mtot.par",
                                parallel="Mrad.par")
  
  levels(m$variable) <- list(Mtot="Mtot.perp",
                             Mtot="Mtot.par",
                             Mrad="Mrad.perp",
                             Mrad="Mrad.par")
  invisible(m)
}

params <- expand.grid(d=c(1,5,10), material=c("silver", "gold"), stringsAsFactors = FALSE)
all <- mdply(params, distance)


p <- 
ggplot(all, aes(wavelength, value, colour=factor(d), linetype=orientation))+
  facet_grid(variable~material, scales="free_y") + 
  geom_path() + labs(colour="distance /nm", y="EM enhancement factor", x="wavelength /nm")+
  scale_y_log10() +
  theme_bw() 

## p
p
