library(planar)
require(ggplot2)
require(plyr)

lambda <- 633
k0 <- 2*pi/lambda
kx <- k0*sin(pi/4)*1.5
epsilon <- c(1.5^2, epsAg(lambda)$epsilon, 0.093284^2, 1.0^2, 1.0^2)
thickness <- c(0, 50,1, 10, 0)
angle <- 0.7617526
# angle <- 0
# angle <- 5*pi/180
psi <- 0
psi <- pi/2
psi <- 0
psi <- pi/3

struct <- tamm_stack(wavelength = 704,
                     lambda0=650, N=6, incidence = "left",
                     nH = 1.7, nL = 1.3, dm = 28, 
                     nleft = 1.5, nright = 1.0)
epsilon <- epsilon_dispersion(struct$epsilon, 700)

# a <- cpp_multilayer_field(2*pi/700, 0, unlist(epsilon), 
#                         struct$thickness, seq(-100, 1000), psi)
# 
# b <- multilayerfull(wavelength=700, angle=angle, 
#                        epsilon=epsilon,
#                        thickness=struct$thickness, psi=psi, z=seq(-100, 1000))
# 
# a$ts
# b$ts
# 
# a$tp
# b$tp
res <- internal_field(wavelength = 700, 
               angle = angle, 
               psi = psi, 
               thickness = struct$thickness, 
               dmax = 100, res = 100,  res2=1e4,
               epsilon = unlist(epsilon_dispersion(struct$epsilon, 700)))
              

test <- multilayerfull(wavelength=700, angle=angle, 
                       epsilon=epsilon_dispersion(struct$epsilon, 700),
                       thickness=struct$thickness, psi=psi, z=seq(-100, 1000))

testing <- data.frame(x=seq(-100, 1000), I=test$I)

limits <- ddply(res, .(id), summarize,  material=unique(material),
                xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-Inf, 
                ymax=Inf)


ggplot(res) +
  geom_rect(data=limits, alpha=0.2, colour=NA,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                fill=material))+
  geom_line(aes(x, I)) +
  geom_line(aes(x, I), data=testing, col="red") +
  theme_minimal() + theme(panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank()) +
  guides(fill = guide_legend(override.aes=list(size=1)))


