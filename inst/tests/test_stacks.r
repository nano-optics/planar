library(planar)
a <- dbr_stack()
b <- tamm_stack()
par(mfrow=c(1,2))
plot(a)
plot(rev(a))

require(ggplot2)
autoplot(a)
autoplot(c(a, rev(a)))

autoplot(b)

ff <- simulate_ff(b, wavelength=seq(600, 1000))

ggplot(ff, aes(wavelength, A)) +
  geom_line()

optimum <- subset(ff, A == max(A))
nf <- simulate_nf(b, wavelength=optimum$wavelength)

autoplot(b) +
  geom_line(aes(x, I), data=nf)
