
lf <- list.files(pattern = "rmd$")
#plyr::l_ply(lf, knitr::knit2html)
plyr::l_ply(lf, knitr::purl)
#plyr::l_ply(gsub("rmd$", "r", lf), file.copy, "../../demo/", overwrite=TRUE)

