##' radiationGrob
##'
##' radiationGrob
##' @title radiationGrob
##' @export
##' @aliases grid.radiation radiationGrob
##' @param theta angle
##' @param r radius
##' @param start starting angle
##' @param log logical, log-transformation of radius
##' @return grob 
##' @author baptiste Auguie
radiationGrob <- function(theta, r, start=pi, log=TRUE){

  if(log){
    r <- log10(r)
    r <- r - min(r)
  }
  x <- r*cos(theta+start)
  y <- r*sin(theta+start)
  
  g <- linesGrob(x,y, default.units = "native")
  g2 <- linesGrob(range(x), 0*range(y), gp=gpar(col="red"),
                  default.units = "native")
  gTree(children=gList(g, g2), vp=dataViewport(xD=x, yD=y))
}

grid.radiation <- function(...)
  grid.draw(radiationGrob(...))
