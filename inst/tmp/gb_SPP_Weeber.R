library(planar)
library(plyr)
library(ggplot2)

wavelength <- 800
metal <- (0.180 + 5.12i)^2
epsilon <- list(1.5^2, metal, 1.0)
thickness <- c(0, 50, 0)

## first, check the plane wave result
results <- multilayer(epsilon=epsilon,
                      wavelength=wavelength, thickness=thickness, d=1,
                      angle=seq(0, pi/2, length=2e3), polarisation='p')

maxi <- max(results$Mr.perp[[2]] + results$Mr.par[[2]], na.rm=T)
spp <- results$angle[which.max(results$Mr.perp[[2]] + results$Mr.par[[2]])]
print(spp)

simulation <- function(w0=10){
  w0 <- w0*1e3
  xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0+5000,length=100), y=0, z=thickness[2]+1))
  res <- adply(xyz, 1, gaussian_near_field2, 
               epsilon=unlist(epsilon), thickness=thickness,
               w0=w0, alpha=spp, maxEval=1000)
  data.frame(xyz, field=res[[2]])
}


params <- data.frame(w0=c(10))
all <- mdply(params, simulation, .progress="text")

peak <- subset(all, field == max(field))
peakx <- peak$x /1000
peakl <- round(peakx,2)
peaky <- peak$field

p <- ggplot(all, aes(x/1000, field, group=w0, colour=factor(w0)))+
  geom_line()  +
  geom_vline(aes(x=0,y=NULL),lty=2) +
  annotate("text", label=peakl, y=peaky, x=peakx, hjust=0, vjust=0, fontface="italic") +
  labs(x=expression(x), y=expression("|E|"^2), 
       colour=expression(w[0]/mu*m)) +
  guides(colour="none") +
  theme(panel.background=element_rect(fill=NA)) +
  theme()

print(p)
