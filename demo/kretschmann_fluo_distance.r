## Fluorescence enhancement vs distance in the Kretschmann configuration

library(planar)
library(ggplot2)

## from left to right
## incident glass | metal | water


## variable  plasmon
## 1     gold 54.88244
## 2   silver 53.08154


fluorescence.enhancement <- function(d=seq(1,10),
                                     lambdaex=632.8,
                                     shift=1000,
                                     thetaex=seq(50,55)*pi/180,
                                     thetafluo=thetaex,
                                     nPrism = 1.766,
                                     nWater = 1.33,
                                     material = "silver",
                                     polarisation="p",
                                     thickness = c(0, 50, 0),
                                     Nquadrature1 = 200, Nquadrature2 = 500,
                                     Nquadrature3 = 500,
                                     qcut = NULL, rel.err= 1e-1,
                                     GL=TRUE){

  ## local field enhancement for the excitation of the fluorophore

  lambdafluo <- raman.shift(lambdaex, shift)

  epsilon.ex <- if(material == "silver") epsAg(lambdaex)$epsilon else epsAu(lambdaex)$epsilon 
  epsilon.fluo <- if(material == "silver") epsAg(lambdafluo)$epsilon else epsAu(lambdafluo)$epsilon 

  Mex <- multilayer(lambda = lambdaex, theta = thetaex, polarisation = polarisation,  
                    thickness = thickness, d = d,
                    epsilon = list(nPrism^2, epsilon.ex, nWater^2))
  
  Mfluo <- multilayer(lambda = lambdafluo, theta = thetafluo,
                      polarisation = polarisation,  
                      thickness = thickness, d=d,
                      epsilon = list(nPrism^2, epsilon.fluo, nWater^2))

  if(GL){
    
    params <- list(lambda = lambdafluo,
                   epsilon = list(nWater^2, epsilon.fluo, nPrism^2), # reversed order
                   thickness = thickness,
                   Nquadrature1 = Nquadrature1, Nquadrature2 = Nquadrature2,
                   Nquadrature3 = Nquadrature3, qcut = NULL)
    
  Mtot <-   sapply(d, function(.d){
    dl <- do.call(dipole.GL, c(params, list(d=.d)))
    dl
  }, simplify=TRUE)
    
  } else {
    params <- list(lambda = lambdafluo,
                   epsilon = list(nWater^2, epsilon.fluo, nPrism^2), # reversed order
                   thickness = thickness,
                   Nquadrature1 = Nquadrature1, Nquadrature2 = Nquadrature2,
                   Nquadrature3 = Nquadrature3, qcut = qcut, rel.err=rel.err)
    
    Mtot <-   sapply(d, function(.d){
      dl <- do.call(dipole, c(params, list(d=.d, show.messages = TRUE)))
      ## print(dl$evaluations)
      dl$results
    }, simplify=TRUE)
    
  }
  
  
  ## incident left -> glass | metal | water
  last <- length(Mex$Mr.perp)
  Mex.perp <-  Mex$Mr.perp[[last]]
  Mfluo.perp <-  Mfluo$Mr.perp[[last]]
  Mex.par <-  Mex$Mr.par[[last]]
  Mfluo.par <-  Mfluo$Mr.par[[last]]
  
  Mtot.perp <- matrix(unlist(Mtot["Mtot.perp",]),
                      nrow=nrow(Mex.perp), ncol=ncol(Mex.perp), byrow=TRUE)
  Mtot.par <- matrix(unlist(Mtot["Mtot.par",]),
                     nrow=nrow(Mex.perp), ncol=ncol(Mex.perp), byrow=TRUE)
  
  d <- data.frame(d=rep(d, each=length(thetaex)),
                  thetaex=rep(thetaex, length(d)),
                  Mex.perp = c(Mex.perp),
                  Mex.par = c(Mex.par),
                  Mfluo.perp = c(Mfluo.perp),
                  Mfluo.par =  c(Mfluo.par),
                  Mtot.perp = c(Mtot.perp),
                  Mtot.par = c(Mtot.par),
                  F.perp = c(Mex.perp * Mfluo.perp / Mtot.perp),
                  F.par = c(Mex.par * Mfluo.par / Mtot.par))

  
  m <- melt(d, id=c("d", "thetaex"))
  ## reorder the levels
  levels(m$variable) <- c("Mex.perp", "Mex.par", "Mfluo.perp", "Mfluo.par",
                          "Mtot.perp", "Mtot.par", "F.perp", "F.par")
  m$type <- factor(m$variable)
  levels(m$type) <- list("Mex"="Mex.perp", "Mfluo"="Mfluo.perp", "Mtot"="Mtot.perp",
                         "Mex"="Mex.par", "Mfluo"="Mfluo.par", "Mtot"="Mtot.par",
                         "F"="F.perp", "F"="F.par")
  m$orientation <- factor(m$variable)
  levels(m$orientation) <- list("perp"="Mex.perp", "perp"="Mfluo.perp", "perp"="Mtot.perp",
                                "par"="Mex.par", "par"="Mfluo.par", "par"="Mtot.par",
                                "perp"="F.perp", "par"="F.par")
  
  m
}


system.time(
test <- fluorescence.enhancement(d=seq(1,50, by=1), thetaex=seq(51,53, by=0.1)*pi/180,
                                 Nquadrature1 = 15, Nquadrature2 = 135,
                                 Nquadrature3 = 135)
            )

## system.time(
## test2 <- fluorescence.enhancement(d=seq(1,50, by=1), thetaex=c(53,53)*pi/180, GL=FALSE,
##                                   Nquadrature1 = 15, Nquadrature2 = 135,
##                                  Nquadrature3 = 135)
##             )


m2 <- subset(test, !variable %in% c("Mfluo.par", "Mfluo.perp"))
m2$value[m2$type == "Mtot" ] <- log10(m2$value[m2$type == "Mtot" ])

modify.levels <- function(f, modify=list()){
  f <- factor(f)
 levs = levels(f)
 m = match(modify,levs)
 levs[m] = names(modify)
 factor(f,labels=levs)
}

m2$variable <- modify.levels(m2$variable,
                             list("log.Mtot.perp"="Mtot.perp",
                                  "log.Mtot.par"="Mtot.par"))
p <- 
qplot(d, value, data=m2, color=thetaex*180/pi, group=thetaex*180/pi, geom="line") +
  facet_wrap(~variable,scales="free", ncol=2,as.table=FALSE)+
  scale_x_continuous("distance /nm", expand=c(0,0), lim=c(0, max(test$d)))+
  scale_y_continuous("")+ labs(colour=expression(theta[laser]))+
  geom_hline(yintercept=0)+
  theme_minimal()

p
