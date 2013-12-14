##' Multilayer optics
##'
##' R/c++ implementation of the dipole emission near a planar multilayer stack
##' 
##' @name planar-package
##' @docType package
##' @useDynLib planar
##' @import Rcpp
##' @import RcppArmadillo
##' @import dielectric
##' @import methods
##' @import plyr 
##' @import reshape2 
##' @importFrom statmod gauss.quad
##' @importFrom cubature adaptIntegrate
##' @title planar
##' @author baptiste Auguie \email{baptiste.auguie@@gmail.com}
##' @references
##' Etchegoin, P. Le Ru, E., Principles of Surface-Enhanced Raman Spectroscopy, Elsevier, Amsterdam (2009).
##' 
##' L. Novotny, E. Hecht, Principles of Nano-optics Cambridge University Press, 2006
##' 
##' H. Raether. Surface Plasmons on Smooth and Rough Surfaces and on Gratings. Springer, 1988.
##' @keywords packagelibrary
##' 
NULL

##' Rcpp module: planar
##' 
##' Exposes C++ functions multilayer and recursive_fresnel
##' @name planar
##' @docType data
##' @export
##' @details
##' \itemize{
##'  \item{multilayer}{reflection and transmission of a multilayer using a transfer matrix formalism}
##'  \item{recursive_fresnel}{reflection and transmission of a multilayer using recursive application of Fresnel coefficients}
##' }
##' @examples
##' show( planar )
NULL

##' Rcpp module: gaussian
##' 
##' Exposes C++ function integrand_gb
##' @name gaussian
##' @docType data
##' @export
##' @details
##' \itemize{
##'  \item{integrand_gb}{integrand for gaussian beam excitation at a planar interface}
##' }
##' @examples
##' show( gaussian )
NULL