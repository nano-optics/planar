##' minimal theme
##'
##' bare-bones theme for ggplot2
##' @title theme_minimal
##' @export
##' @param base_size font size
##' @param base_family font family
##' @return theme (list)
##' @author baptiste Auguie
##' 
theme_minimal <- 
function (base_size = 12, base_family = "") 
{
    structure(list(axis.line = theme_blank(), axis.text.x = theme_text(family = base_family, 
        size = base_size * 0.8, lineheight = 0.9, vjust = 1), 
        axis.text.y = theme_text(family = base_family, size = base_size * 
            0.8, lineheight = 0.9, hjust = 1), axis.ticks = theme_segment(colour = "black", 
            size = 0.2), axis.title.x = theme_text(family = base_family, 
            size = base_size, vjust = 1), axis.title.y = theme_text(family = base_family, 
            size = base_size, angle = 90, vjust = 0.5), axis.ticks.length = unit(0.3, 
            "lines"), axis.ticks.margin = unit(0.5, "lines"), 
        legend.background = theme_rect(colour = NA),
                   legend.key = theme_rect(colour = NA), 
        legend.key.size = unit(1.2, "lines"), legend.key.height = NA, 
        legend.key.width = NA, legend.text = theme_text(family = base_family, 
            size = base_size * 0.8), legend.text.align = NA, 
        legend.title = theme_text(family = base_family, size = base_size * 
            0.8, face = "bold", hjust = 0), legend.title.align = NA, 
        legend.position = "right", legend.direction = "vertical", 
        legend.box = NA, panel.background = theme_rect(fill = "white", 
            colour = NA), panel.border = theme_rect(fill = NA, 
            colour = "grey90"), panel.grid.major = theme_line(colour = "grey90", 
            size = 0.2), panel.grid.minor = theme_line(colour = "grey98", 
            size = 0.5), panel.margin = unit(0.25, "lines"), 
        strip.background = theme_rect(fill = NA, colour = NA), 
        strip.text.x = theme_text(family = base_family, size = base_size * 
            0.8), strip.text.y = theme_text(family = base_family, 
            size = base_size * 0.8, angle = -90), plot.background = theme_rect(colour = NA), 
        plot.title = theme_text(family = base_family, size = base_size * 
            1.2), plot.margin = unit(c(1, 1, 0.5, 0.5), "lines")), 
        class = "options")
}
