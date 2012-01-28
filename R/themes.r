##' minimal theme
##'
##' bare-bones theme for ggplot2
##' @title theme_minimal
##' @param base_size font size
##' @param base_family font family
##' @return theme (list)
##' @author baptiste Auguié
##' 
theme_minimal <- 
function (base_size = 12, base_family = "") 
{
    structure(list(axis.line = theme_blank(), axis.text.x = theme_text(family = base_family, 
        size = base_size * 0.8, lineheight = 0.9, vjust = 1), 
        axis.text.y = theme_text(family = base_family, size = base_size * 
            0.8, lineheight = 0.9, hjust = 1), axis.ticks = theme_segment(colour = "black", 
            size = 0.2), axis.title.x = theme_text(family = base_family, 
            size = base_size, vjust = 0.5), axis.title.y = theme_text(family = base_family, 
            size = base_size, angle = 90, vjust = 0.5), axis.ticks.length = unit(0.15, 
            "cm"), axis.ticks.margin = unit(0.1, "cm"), legend.background = theme_rect(colour = NA, fill=NA), 
        legend.margin = unit(0.2, "cm"), legend.key = theme_blank(), 
        legend.key.size = unit(1.2, "lines"), legend.key.height = NULL, 
        legend.key.width = NULL, legend.text = theme_text(family = base_family, 
            size = base_size * 0.8), legend.text.align = NULL, 
        legend.title = theme_text(family = base_family, size = base_size * 
            0.8, face = "bold", hjust = 0), legend.title.align = NULL, 
        legend.position = "right", legend.direction = NULL, legend.justification = "center", 
        legend.box = NULL, panel.background = theme_blank(), panel.border = theme_blank(),
                   panel.grid.major = theme_line(colour = "grey90", 
            size = 0.2), panel.grid.minor = theme_line(colour = "grey98", 
            size = 0.5), panel.margin = unit(0.25, "lines"), 
        strip.background = theme_blank(), 
        strip.text.x = theme_text(family = base_family, size = base_size * 
            0.8), strip.text.y = theme_text(family = base_family, 
            size = base_size * 0.8, angle = -90), plot.background = theme_blank(), 
        plot.title = theme_text(family = base_family, size = base_size * 
            1.2), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")), 
        class = "options")
}
