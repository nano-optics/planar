## Integrated decay rates and efficiency for a dipole near a semi-infinite air/metal interface
## for gold and silver, varying the wavelength and the dipole-interface distance

library(planar)
library(ggplot2)
library(lattice)
require(reshape2)
library(gridExtra)
require(plyr)

wvl <- seq(250, 750, by=250) 

silver <- epsAg(wvl)
gold <- epsAu(wvl)


distance <- function(d, material="silver", ...){

  material <- get(material)
  
  params <- list(d=d,
                 lambda = material$wavelength,
                 epsilon = list(incident=1.0^2, material$epsilon),
                 thickness = c(0, 0))
  
  res <- do.call(dipole, c(params, list(Nquadrature1 = 0, Nquadrature2 = 0,
                                        Nquadrature3 = 0, qcut = NULL,
                                       rel.err=1e-2, show.messages=FALSE)))

  res$Q.perp <- res$Mrad.perp / res$Mtot.perp
  res$Q.par <- res$Mrad.par / res$Mtot.par
  
  m <- melt(res, id = "wavelength")
  ## print(paste(c("evaluations", dl$evaluations), collapse=" "))
  ## print(paste(c("Max rel. error", round(100*max(unlist(dl$errors)),1), "%"), collapse=" "))
  
  m$orientation <- m$variable
  
  levels(m$orientation) <- list(perpendicular="Mtot.perp",
                                perpendicular="Mrad.perp",
                                perpendicular="Q.perp",
                                parallel="Mtot.par",
                                parallel="Mrad.par",                                
                                parallel="Q.par")
  
  levels(m$variable) <- list(Mtot="Mtot.perp",
                             Mtot="Mtot.par",
                             Mrad="Mrad.perp",
                             Mrad="Mrad.par",
                             Q="Q.perp",
                             Q="Q.par")
  invisible(m)
}

params <- expand.grid(d=seq(2, 50, by=2),
                      material=c("silver", "gold"), stringsAsFactors = FALSE)
progress <- if(interactive()) "text" else "none"
all <- mdply(params, distance, .progress=progress)


p <- 
  ggplot(all, aes(d, value, linetype=orientation, color=variable))+
  facet_grid(material~wavelength, scales="free_y") + 
  geom_line() + labs(colour="variable", y="EM enhancement factor", x="d /nm")+
  scale_y_log10() +
  theme(legend.direction="horizontal", legend.position="bottom")

p
