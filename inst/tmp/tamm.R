require(planar)
require(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA),
          panel.margin=unit(4,"mm")))

.col <- c("black", RColorBrewer::brewer.pal(6,"Set2"))
.fill <- c("white", RColorBrewer::brewer.pal(6,"Pastel2"))
col_palette <- .col[c(4,7,3,5,6,1)]
fill_palette <- .fill[c(4,7,3,5,6,1)]
# grid.raster(.fill,int=F)
          
## ----bottomstructure----
pars <- list(lambda0=630, 
             angle=0, incidence="right",
             N=8, dx1=0, dm=20,
             n1=1.28, n2=1.72, dbleft=0,
             nleft=1.0, nright=1.5)

dL <- pars$lambda0/4/pars$n1
dH <- pars$lambda0/4/pars$n2

p.nf <- within(pars, { wavelength <- 700; 
                       res <- 1e4;
                       dmax <- 200;
                       fun  <-  tamm_stack})

s <- do.call(tamm_stack, p.nf)

d.nf <- do.call(planar:::simulate_nf, p.nf)

d.nf <- recode_dbr(d.nf, before=c("substrate","", "Au",""), 
                   after=c("", "superstrate"), dbr = c("nH", "nL"))

limits <- ddply(d.nf, .(id), summarize, abb=unique(abb),
                material=unique(material),
                xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-1, 
                ymax=1)

  ggplot(d.nf) + 
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax,
                group=id), data=limits, fill="white") +
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax-0.2,
                fill=material, group=id), col=NA,  alpha=0.5,data=limits) +
  geom_text(aes(x=xmid, y=ymax, label=abb, group=id,colour=material), size=4,
            hjust=0.5, vjust=1,angle=0,
            data=limits)+
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous("", breaks=NULL) +
  scale_colour_manual("", values=col_palette) +
  scale_fill_manual("", values=fill_palette) +
  guides(colour="none",fill="none",alpha="none")+
  theme_minimal()


## ----dx1ff----
params <- expand.grid(dm=c(0, 20),
                      dx1 = seq(-40, 40, by=20), 
                      incidence=c("left", "right"),
                      stringsAsFactors=FALSE)

d1 <- mdply(params, ff_simulation, fun=porous_stack, 
            lambda0=630, 
            angle=0, 
            N=8, #dx1=0, dm=20,incidence="right",
            n1=1.28, n2=1.72, dbleft=0,
            nleft=1.0, nright=1.5)

m1 <- melt(d1, id=c("wavelength", "dm", "dx1", "incidence"), 
           meas=c("T", "R", "A"))
m1$lab <- factor(m1$variable, 
                 labels=c("Transmittance", "Reflectance", "Absorbance"))

ggplot(m1, aes(wavelength, value, linetype=factor(dm),
               colour=factor(dx1), group=interaction(dm, dx1))) + 
  geom_line()+ facet_grid(lab~incidence)+
  scale_y_continuous(lim=c(0,1), expand=c(0,0))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_colour_brewer("dx /nm",palette="Set1")+
  labs(x="wavelength /nm", y="")+
  guides(colour="legend", linetype="none")


## ----dx1nf----
minima <- ddply(subset(d1, dm > 0 & wavelength > 550 ), 
                .(dm, dx1, incidence), summarise,
                wavelength=wavelength[R==min(R)], mini=min(R))

d2 <- mdply(minima, nf_simulation, 
              fun=porous_stack, 
            lambda0=630, 
            angle=0, 
            N=8, #dx1=0, dm=20,incidence="right",
            n1=1.28, n2=1.72, dbleft=0,
            nleft=1.0, nright=1.5,
            dbright=200,
            full.thickness=1500)

limits <- ddply(d2, .(id, dm, dx1, incidence), summarize, 
                material=unique(material),
                xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-Inf, 
                ymax=Inf, wavelength=unique(wavelength))

titles <- ddply(limits, .(dm, dx1, incidence), summarise,
                x = 600, y=Inf,
                wavelength=unique(wavelength))

ggplot(d2) + facet_grid(dx1~incidence, labeller=label_both)+
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, 
                group=id), fill="white", data=limits) +
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, 
                fill=material, group=id), data=limits, alpha=0.5) +
  geom_text(data=titles, vjust=1, parse=TRUE,
            aes(x=x,y=y, label=sprintf("lambda == %.2f*nm", wavelength)))+
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous(expression("|E|"^2),expand=c(0,0)) +
  geom_hline(yintercept=0) +
  geom_line(aes(x, I))+
  scale_fill_manual("", values=rev(fill_palette)) +
  guides(colour="none",fill="none")


## ----bottomlown----


objective <- function(p){
  
  ff <- ff_simulation(fun=porous_stack, 
                      dm=p[1], dx2 = p[2], 
                      incidence="right",
                      wavelength = seq(600, 800,length=1e3),
                      n1=1.72, n2=1.28, N=8, lambda0=630,
                      nright=1.5, nleft=1.0)
  min(ff$R)
}

best <- optim(par=c(20, 630/4/1.28), objective,
              lower=c(10, 0), upper=c(50, 200), method="L-BFGS-B")

params <- expand.grid(dm=c(0, best$par[1]),
                      dx2 = best$par[2], 
                      incidence=c("left"),
                      stringsAsFactors=FALSE)

params <- expand.grid(dm=c(0, 20),
                      dx2 = 630/4/1.28, 
                      incidence=c("right"),
                      stringsAsFactors=FALSE)

d1 <- mdply(params, ff_simulation, fun=porous_stack, 
            wavelength = seq(400, 1000,length=1e3),
            n1=1.72, n2=1.28, N=8, lambda0=630,
            nright=1.0, nleft=1.5)

ddply(subset(d1, dm > 0 & wavelength > 550 ), 
      .(dm, dx2, incidence), summarise,
      wavelength=wavelength[R==min(R)], mini=min(R))

m1 <- melt(d1, id=c("wavelength", "dm", "dx2", "incidence"), 
           meas=c("T", "R", "A"))
m1$lab <- factor(m1$variable, 
                 labels=c("Transmittance", "Reflectance", "Absorbance"))

p1 <- ggplot(m1, aes(wavelength, value, linetype=factor(dm),
               colour=factor(dx2), group=interaction(dm, dx2))) + 
  geom_line()+ facet_grid(lab~incidence)+
  scale_y_continuous(lim=c(0,1), expand=c(0,0))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_colour_brewer("dx /nm",palette="Set1")+
  labs(x="wavelength /nm", y="")+
  guides(colour="none", linetype="none")

minima <- ddply(subset(d1, dm > 0 & wavelength > 550 ), 
                .(dm, dx2, incidence), summarise,
                wavelength=wavelength[R==min(R)], mini=min(R))

d2 <- mdply(minima, nf_simulation, 
            fun=porous_stack,  
            n1=1.72, n2=1.28, N=8, lambda0=630,
            nright=1.0, nleft=1.5, dbleft=500)

limits <- ddply(d2, .(id, dm, dx2, incidence), summarize, 
                material=unique(material),
                xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-Inf, 
                ymax=Inf, wavelength=unique(wavelength))

titles <- ddply(limits, .(dm, dx2, incidence), summarise,
                x = 600, y=Inf,
                wavelength=unique(wavelength))

p2 <- ggplot(d2) + #facet_grid(dx2~incidence, labeller=label_both)+
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, 
                group=id), fill="white", data=limits) +
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, 
                fill=material, group=id), data=limits, alpha=0.5) +
  geom_text(data=titles, vjust=1, parse=TRUE,
            aes(x=x,y=y, label=sprintf("lambda == %.2f*nm", wavelength)))+
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous(expression("|E|"^2),expand=c(0,0)) +
  geom_hline(yintercept=0) +
  geom_line(aes(x, I))+
#   scale_fill_manual("", values=rev(fill_palette[-1])) +
  guides(colour="none",fill="none")

grid.arrange(p1,p2, ncol=2)

