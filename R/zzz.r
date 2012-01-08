
planar <- new( "Module" )
 .onLoad <- function(pkgname, libname){
     loadRcppModules(direct=FALSE)
 }

