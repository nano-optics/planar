
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----setup, results='hide'----------------------------------------------------
library(planar)
library(ggplot2)
require(reshape2)
require(plyr)


nf_model <- function(n, k){
  
  
  m <- lfief(wavelength=500, angle=0*pi/180, polarisation='p',
             thickness=c(0, 0), dmax=250, res=5000,
             epsilon=list(1.0^2,  (n+1i*k)^2), 
             displacement=FALSE)
  m$L1 <- factor(m$L1, labels=c("glass",  "air"))
  m
}

field_model <- function(n,k){
  
  res <- internal_field(wavelength=500, angle=0, psi=0,
                        thickness = c(0,0.001,0), 
                        dmax=250,  res=1e3,
                        epsilon=c(1.0^2, 1.0, (n+1i*k)^2), 
                        field = TRUE)
  
  mutate(res, real = Re(Ex), imag = Im(Ex), intensity=Mod(Ex)^2)
}

ff_model <- function(n,k){
  
  res <- recursive_fresnelcpp(epsilon=list(1.0^2, 1.0, (n+1i*k)^2),
                                  wavelength=500, thickness=c(0, 0.001, 0),
                                  angle=0, polarisation='p')

  data.frame(R = res$R, A=1-res$R)
}






params <- expand.grid(n=seq(1,5,length=5), k=c(0.0001,0.1,1,2,10))
m <- mdply(params, nf_model)
m2 <- mdply(params, field_model)
m3 <- mdply(params, ff_model)
m3 <- ddply(m3, c("k"), mutate, y=seq(3.8, 1.5, length=5), label=sprintf("%.3f",A))

## ----plot, fig.width=12--------------------------------------------------

limits <- ddply(m, .(L1), summarize, 
                xmin=min(x), xmax=max(x), ymin=-Inf, ymax=Inf)

ggplot(subset(m, variable=="M.par")) +
  facet_grid(k~., labeller = label_both)+
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, 
                fill=factor(L1)),     data=limits, alpha=0.2) +
  geom_path(aes(x, value, colour=n, group=interaction(L1, n))) +
  geom_text(data=m3, aes(x=120, y=y, label=label, colour=n), hjust=0, show.legend = FALSE) +
  annotate("text", x=100, hjust=1, y=2.5, label="A=1-R=") +
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous(expression("|E|"^2),expand=c(0,0),lim=c(0,4)) +
  geom_hline(yintercept=0) +
  # scale_colour_brewer("", palette="Set1")+
  scale_fill_manual(values=c("white", "grey80")) +
  guides(fill="none", colour=guide_legend()) +
  theme_bw() +
  theme(legend.position="top", panel.margin=unit(1, "line"),
        strip.background=element_blank())

# ggsave("nk.pdf", width=4, height=8)

mm <- melt(m2, id=c("x", "id", "n", "k"), meas=c("real","imag", "intensity"))
ggplot(mm) +
  facet_grid(k~variable, labeller = label_both)+
  geom_line(aes(x, value, colour=n, group=n)) +
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous("E",expand=c(0,0),lim=c(-2,4)) +
  geom_vline(xintercept=0, lty=2) +
  geom_hline(yintercept=0, lty=2) +
  # scale_colour_brewer("", palette="Set1")+
  scale_fill_manual(values=c("white", "grey80")) +
  guides(fill="none", colour=guide_legend()) +
  theme_minimal() +
  theme(legend.position="top")


