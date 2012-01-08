##' transmission loss through a prism
##'
##' transmission loss through a prism
##' @title transmission
##' @export
##' @param nPrism prism refractive index
##' @param extAngle_r external incident angle in radians
##' @param polar polarisation
##' @return transmission 
##' @author baptiste Auguie
transmission <- function(nPrism, extAngle_r, polar = "p"){
  
  alpha <- asin(sin(extAngle_r) / nPrism) # % refracted angle
  
  if(polar == 'p'){
          4 * nPrism*cos(extAngle_r)*cos(alpha) /
            (nPrism * cos(extAngle_r) + cos(alpha))^2
        } else {
         .NotYetImplemented()
        }			
  
}
