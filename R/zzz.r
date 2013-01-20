
planar <- new( "Module" )
gaussian <- new( "Module" )

.onLoad <- function(libname, pkgname) {
     loadRcppModules(direct=FALSE)
}
