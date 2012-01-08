## radiation pattern of a dipole near a dielectric/(metal)/dielectric interface
## parallel and perpendicular orientations, p- and s- polarisations

library(planar)
library(ggplot2)

## from left to right
## incident glass | (metal) | water

## variable  plasmon
## 1     gold 54.88244
## 2   silver 53.08154

angular.pattern <- function(d=1,
                            lambdaex=632.8,
                            shift=1000,
                            nPrism = 1.766,
                            nWater = 1.33,
                            material = "silver",
                            polarisation="p",
                            thickness = 20){

  lambdafluo <- raman.shift(lambdaex, shift)  
  epsilon.fluo <- if(material == "silver") epsAg(lambdafluo)$epsilon else epsAu(lambdafluo)$epsilon 

  critical.angle <- asin(nWater / nPrism)*180/pi
  
  spp.angle <- asin(Re(sqrt((epsilon.fluo * nWater^2) /
                            (epsilon.fluo + nWater^2)) / nPrism)) *180/pi
  
  epsilon.fluo <- if(is.null(thickness)) NULL else epsilon.fluo
  
  theta <- seq(-90, 90, length=1e3) * pi/180
  
  Mfluo.bottom <- multilayer(lambda = lambdafluo, theta = theta, polarisation = polarisation, 
                             thickness = c(0, thickness, 0), d = d,
                             epsilon = list(nPrism^2, epsilon.fluo, nWater^2))
  
  Mfluo.top <- multilayer(lambda = lambdafluo, theta = theta, polarisation = polarisation,  
                          thickness = c(0, thickness, 0), d = d,
                          epsilon = list(nWater^2, epsilon.fluo, nPrism^2))

  ## look at field in first and last media
  last <- if(is.null(thickness)) 1 else  2 + length(thickness) - 1
  theta=c(pi+asin(Mfluo.top$q), asin(Mfluo.bottom$q))*180/pi+90
  intensity.par = c(Mfluo.top$Ml.par[[1]] , Mfluo.bottom$Mr.par[[last]])
  intensity.perp = c(Mfluo.top$Ml.perp[[1]] , Mfluo.bottom$Mr.perp[[last]])
  id=rep(1:2, each=length(Mfluo.bottom$R))
  
  res <- 
    data.frame(theta=theta, parallel=intensity.par,
               perpendicular=intensity.perp, id=id )
  
  list(results=melt(res[complete.cases(res),], id=c("theta", "id")),
       critical.angle=critical.angle,
       spp.angle=spp.angle)
}


plot.pattern <- function(m, spp=TRUE, log=TRUE, ...){

  my.labels <- c(45, 0, 45, 90, 45, 0, 45, 90,-90)
  
  ## m$results$value[m$results$value <= 2] <- 2
  ## print(range(m$results$value))
  spp <- if(spp)   geom_path(aes(x=x,y=y, group=id),
                             data=data.frame(x=rep(90 + m$spp.angle *c(-1,1), each=2),
                               id=rep(1:2, each=2),
                               y=rep(range(m$results$value)), 2),
                             colour="grey", linetype=3) else NULL

  pp <- 
    ggplot(m$results) + coord_polar(start=pi/2)+ xlim(0,360) +
      spp +
        geom_path(aes(x=x,y=y, group=id),
                  data=data.frame(x=rep(90 + m$critical.angle *c(-1,1), each=2),
                    id=rep(1:2, each=2),
                    y=rep(range(m$results$value)), 2),
                colour="grey", linetype=2) +  
                  geom_path(aes(x=x,y=y, group=id),
                            data=data.frame(x=c(180, 180, 360, 360), id=rep(1:2, each=2),
                              y=rep(range(m$results$value), 2)),
                            colour="grey", linetype=1, size=2) +  
                              geom_polygon(aes(theta, value, group=interaction(id, variable), fill=factor(id)), alpha=0.1) +
                                  geom_path(aes(theta, value, group=interaction(id, variable),
                                                colour=factor(id), linetype=variable)) +
                                                  scale_x_continuous(breaks=seq(0,360,by=45), minor_breaks=seq(0,360,by=15), labels=my.labels, expand=c(0,0))+
  scale_fill_brewer(palette="Pastel1", legend=FALSE)+
  scale_colour_brewer(palette="Set1", legend=FALSE)+
  theme_minimal() +
  labs(x="", y="", linetype="orientation") +
  opts(legend.position="top", legend.direction="horizontal", ...)

 if(log) pp + scale_y_log10() else pp
}

library(gridExtra)

metalp <- angular.pattern(lambdaex=632.8, pol="p", thick=50, d=10)
dielectricp <- angular.pattern(lambdaex=632.8, pol="p", thick=0, d=10)

metals <- angular.pattern(lambdaex=632.8, pol="s", thick=50, d=10)
dielectrics <- angular.pattern(lambdaex=632.8, pol="s", thick=0, d=10)



p1 <- plot.pattern(metalp,log=T, title="p-polarisation, metal")
p2 <- plot.pattern(dielectricp, spp=FALSE, log=T, title="p-polarisation, dielectric")
p3 <- plot.pattern(metals,log=F, title="s-polarisation, metal")
p4 <- plot.pattern(dielectrics, spp=FALSE, log=F, title="s-polarisation, dielectric")

if(interactive())
  grid.arrange(p1, p2, p3, p4, ncol=2)
