## radiation pattern of a dipole near a dielectric/(metal)/dielectric interface
## parallel and perpendicular orientations, p- polarisation

library(planar)
require(reshape2)
library(ggplot2)
library(gridExtra)

library(dielectric)

data(AgPalik)
AgPalik$set_span(0.2, 0.9)

wv <- c(seq(200, 300, by=5), seq(305, 400, by=5), seq(410, 900, by=10))

silver <- AgPalik$predict(new.wavelength = wv * 1e-3, all.knots=TRUE)
d <- dielectric2plot(silver)

p2 <- 
 ggplot(d, aes(1e3*wavelength, value)) + geom_path() +
   facet_grid(variable~., scales="free") +
  geom_hline(aes(yintercept=y), lty=2,
             data=data.frame(variable = factor("real", levels = c("real", "imag")), y = -1)) +
  labs(x = expression("Wavelength / nm"), y="Dielectric function")


## wrap simulation in a function to loop over wavelengths
simulation <- function(lambda = 632.8){
  
  eps <- AgPalik$predict(new.wavelength = lambda * 1e-3, all.knots=TRUE)
  current <- dielectric2plot(eps)


  params <- list(lambda=lambda,
                 nPrism = 1.5,
                 nWater = 1.0,
                 metal = eps$epsilon,
                 theta = seq(-90, 90, length=1000) * pi/180)
  
  back <- list(lambda = params$lambda,
               theta = params$theta, 
               epsilon = list(params$nPrism^2, params$metal, params$nWater^2),
               thickness = c(0, 50, 0),
               d = 1,
               polarisation = 'p')
  
  critical <- 270 +c(-1, 1) * asin(1/1.5) *180/pi 
  
  front <- invert_incidence(back)
  
  ## actual simulation
  M.front <- do.call(multilayer, front)
  M.back <- do.call(multilayer, back)
  
  theta <- c(params$theta + pi/2, params$theta + 3*pi/2) * 180 / pi
  ## look at field in first and last media
  last <- length(M.back$Mr.par)
  intensity.par = c(M.front$Ml.par[[1]] , M.back$Mr.par[[last]])
  intensity.perp = c(M.front$Ml.perp[[1]] , M.back$Mr.perp[[last]])
  
  combined <- data.frame(theta=theta, parallel=intensity.par,
                       perpendicular=intensity.perp,
                         side = gl(2, length(M.front$q)))

  
  ## polar plot with two variables
  m <- melt(combined, id=c("theta", "side"))
  lambda <<- lambda 
  p <- 
    ggplot(m, aes(theta, value)) +
      annotate("segment", x=critical[1], xend=critical[1], y=0, yend=1, colour="green", size=0.5) +
        annotate("segment", x=critical[2], xend=critical[2], y=0, yend=1, colour="green", size=0.5) +
          geom_polygon(aes(fill=side), alpha=0.5)  +
            geom_path(aes(linetype=variable), alpha=0.5)  +
              scale_x_continuous(breaks = seq(0, 360, by = 45), limits = c(0, 360), expand = c(0, 0)) +
                scale_y_log10(lim=c(1e-4, 200)) +
                  coord_polar(start  =3*pi/2) +
                    annotate("segment", x=360, xend=360, y=0, yend=1, colour="orange", size=2) +
                      annotate("segment", x=180, xend=180, y=0, yend=1, colour="orange", size=2) +
                        scale_fill_brewer(palette = "Pastel1") +
                          labs(x = "Angle / degrees", y = "LFIEF", fill = "side", linetype = "dipole") +
                            opts(title = bquote(lambda == .(sprintf("%.1f nm", lambda))))
  
  
  perm <- p2 + geom_point(data=current, colour="red", size=2) +
    opts(title = bquote(lambda == .(sprintf("%.1f nm", lambda))))
  
  grid.arrange(p, perm, ncol=2)
}

simulation()


library(animation)

## saveSWF({

##    l_ply(wv, simulation, .progress = "text")
## }, ani.dev = "pdf", ani.type = "pdf", swf.name = "kretschmanns.swf", interval = 0.5, nmax = length(wv), 
##     ani.height = 4, ani.width = 12, outdir = "/Users/auguieba/Dropbox/Public/", autobrowse=FALSE, verbose=FALSE)

## saveGIF({
##        l_ply(wv, simulation, .progress = "text")
## }, movie.name = "kretschmann.gif", interval = 0.1, nmax = 30, ani.width = 1200, 
##     ani.height = 400)
