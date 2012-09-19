##' minimal theme
##'
##' bare-bones theme for ggplot2
##' @title theme_minimal
##' @param base_size font size
##' @param base_family font family
##' @return theme (list)
##' @author baptiste Auguie
##' @export
##' @family user_level theme
theme_minimal <- function(base_size = 12, base_family = "Helvetica") {
  # Starts with theme_bw and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          plot.background = element_blank()
          )
}

if (getRversion() >= "2.15.1") 
  utils::globalVariables(c("%+replace%", 
                           "theme", 
                           "element_blank"))
