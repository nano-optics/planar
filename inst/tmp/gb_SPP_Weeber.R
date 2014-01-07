library(planar)
library(plyr)
library(ggplot2)

simulation <- function(d, probe=50, w0=10e3) {
  
  s <- list(wavelength=800, thickness=c(0,d,0),
                           epsilon=list(1.52^2, (0.180 + 5.12i)^2, 1.0^2))

  pw <- multilayer(epsilon=s$epsilon,
                   wavelength=s$wavelength, thickness=s$thickness, d=probe,
                   angle=seq(0, pi/2, length=2e3), polarisation='p')
  
  maxi <- max(pw$Mr.perp[[2]] + pw$Mr.par[[2]], na.rm=T)
  spp <- pw$angle[which.max(pw$Mr.perp[[2]] + pw$Mr.par[[2]])]
  
  xyz <- as.matrix(expand.grid(x=seq(-50e3, 150e3,length=100), y=0, z=s$thickness[2]+probe))
  
  res <- gaussian_near_field_ml(xyz, wavelength=s$wavelength,
                                epsilon=unlist(s$epsilon), thickness=s$thickness,
                                w0=w0, alpha=spp, tol=1e-3, maxEval=0)
  data.frame(xyz, field=res/max(res))
}

params <- data.frame(d=c(50, 100))
all <- mdply(params, simulation)

bare <- simulation(0, 50)

peak <- function(d){
peak <- subset(d, field == max(field))
peakx <- peak$x /1000
peakl <- round(peakx,2)
peaky <- peak$field
data.frame(peakx=peakx, peakl=peakl, peaky=peaky)
}

ann <- ddply(all, "d", peak)

p <- ggplot(all, aes(x/1000, field))+ facet_grid(d~., scales="free")+
  geom_line()  +
  geom_line(data=bare, linetype="dotted") +
  geom_vline(aes(x=0,y=NULL),lty=2) +
  geom_blank(data=ann, aes(y=peaky*1.1, x=peakx)) +
  geom_text(data=ann, aes(label=peakl, y=peaky, x=peakx), hjust=0, vjust=0, fontface="italic") +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression(x), y=expression("|E|"^2), 
       colour=expression(w[0]/mu*m)) +
  guides(colour="none") +
  theme_minimal()+
  theme(panel.background=element_rect(fill=NA)) +
  theme()

print(p)
