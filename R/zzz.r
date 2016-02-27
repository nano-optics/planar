
planar <- new( "Module" )
gaussian <- new( "Module" )
collection <- new( "Module" )


.onLoad <- function(libname, pkgname) {
     loadRcppModules(direct=FALSE)
}

#' fill_palette
#' 
#' palette
#' @name fill_palette
#' @rdname palette
#' @docType data
#' 
NULL
