library(R.matlab)
setwd("~/Documents/github/tamm/inst/tests")
tmp <- readMat("ff.mat")
ffm <- data.frame(wavelength = c(tmp$wavelength)*1000,
                 front = c(tmp$R0.TM.TM), 
                 back = c(tmp$R0sub.TM.TM),
                 variable="R")

tmp <- readMat("nf1.mat")
nf1m <- data.frame(x = c(tmp$z*1000),
                 I = c(tmp$Intensity[,1]))
tmp <- readMat("nf2.mat")
nf2m <- data.frame(x = c(tmp$z*1000),
                  I = c(tmp$Intensity[,1]))
# plot(nf)

library(tamm)
# prettier palette
palette(palette_tamm)

s1 <- tamm_stack_ir(pairs=15, dm=20, lambda0 = 900, n1=3.0, n2=3.5, nleft = 3.5)
s2 <- tamm_stack_ir(pairs=15, dm=20, incidence="right",
                    lambda0 = 900, n1=3.0, n2=3.5, nleft = 3.5)
autoplot(s2)

(p <- autoplot(s1))
# 
# s1 <- tamm_stack_ir(900, n1 = 3.5, n2=3, pairs = 15, dm = 20, nleft=1.0, nright=3.5,
#                    position = "before", incidence="left")


# ---- ff -----
ff1 <- simulate_ff(s=s1, wavelength=seq(800, 1100))
ff2 <- simulate_ff(s=s2, wavelength=seq(800, 1100))
head(ff1)
mff1 <- melt(ff1, meas=c("R","T","A"))
mff2 <- melt(ff2, meas=c("R","T","A"))

ggplot(mff, aes(wavelength, value, colour=variable))+
  geom_line() +
  geom_line(data=mff2,lty=2) +
  geom_point(data=ffm,aes(wavelength, front)) +
  geom_point(data=ffm,aes(wavelength, back)) +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="")


# ---- nf -----
nf1 <- simulate_nf(s=s1, wavelength=950)
nf2 <- simulate_nf(s=s2, wavelength=950)
head(nf)

ggplot(nf1, aes(x,I)) + geom_line() +
  geom_line(aes(x+150,I),data=nf2m,lty=2) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0)) 

ggplot(nf2, aes(x,I)) + geom_line() +
  geom_line(aes(max(nf2$x)-150-x,I),data=nf1m,lty=2) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))

