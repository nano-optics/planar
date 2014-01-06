
lf <- list.files(pattern = "rmd$")
try_run <- function(x)
  try(knitr::knit2html(x))
plyr::l_ply(lf, try_run)
# plyr::l_ply(lf, knitr::purl)
#plyr::l_ply(gsub("rmd$", "r", lf), file.copy, "../../demo/", overwrite=TRUE)
 
# cat(gsub("\\.rmd$", "      \n", lf), sep="")

