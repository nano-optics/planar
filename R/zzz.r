
planar <- new( "Module" )
gaussian <- new( "Module" )
collection <- new( "Module" )


.onLoad <- function(libname, pkgname) {
     loadRcppModules(direct=FALSE)
}
