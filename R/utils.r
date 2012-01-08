
## Wide to long format data.frame with new factor variable(s) describing the original columns
classify <- function(d, id=NULL, vars=NULL, ...){

  m <- melt(d, id.vars=id, ...)

  id.variables <- list()
  for (ii in seq_along(vars)){
    id.variables[[ii]] <- rep(vars[[ii]], each=nrow(d))
  }
  names(id.variables) <- names(vars)

  data.frame(m, id.variables)
}

Curry <- function (FUN, ...) 
{
    .orig = list(...)
    function(...) do.call(FUN, c(.orig, list(...)))
}
