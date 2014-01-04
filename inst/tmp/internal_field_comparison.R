library(planar)

lambda <- 633
k0 <- 2*pi/lambda
kx <- k0*sin(pi/4)*1.5
epsilon <- c(1.5^2, epsAg(lambda)$epsilon, 0.093284^2, 1.0^2, 1.0^2)
thickness <- c(0, 50,1, 10, 0)
psi <- 0
angle <- 0.7617526
angle <- 5*pi/180
require(tamm)
struct <- tamm_stack(wavelength = 704,
                     lambda0=650, N=6, incidence = "left",
                     nH = 1.7, nL = 1.3, dm = 28, 
                     nleft = 1.5, nright = 1.0)
# 
# struct <- structure(list(epsilon = list(1, -16.1322843565113+1.1013875606378i, 1), 
#                          wavelength = 697, thickness = c(0, 50, 0), angle = 0, polarisation = "p"), .Names = c("epsilon","wavelength", "thickness", "angle", "polarisation"))

comp <- field_profile(wavelength = struct$wavelength, 
               angle = angle, 
               polarisation = "p", 
               thickness = struct$thickness, 
               dmax = 100, res = 100,  res2=1e3,
               epsilon = struct$epsilon, displacement = FALSE)
              

d1 <- subset(comp, variable == "M.perp")
d2 <- subset(comp, variable == "M.par")
d1$value <- d1$value + d2$value

res <- internal_field(wavelength=struct$wavelength, angle=angle, psi=psi,
                           thickness = struct$thickness, 
                           dmax=100,  res=1e3,
                           epsilon=unlist(struct$epsilon), 
                           field = FALSE)


limits <- ddply(d1, .(L1), summarize,   material=unique(material),
                xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-Inf, 
                ymax=Inf)

limits2 <- ddply(res, .(id), summarize,  material=unique(material),
                xmin=min(x), xmax=max(x), xmid=mean(x), ymin=-Inf, 
                ymax=Inf)


ggplot(d1) +
#   geom_rect(data=limits, fill=NA, linetype="dotted", colour="grey50",
#             aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax))+
  geom_rect(data=limits, alpha=0.2, colour=NA,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                fill=material))+
  geom_line(aes(x, value)) +
  geom_line(data=res, aes(x, I), colour="red", linetype="dashed") +
  theme_minimal() + theme(panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank()) +
  guides(fill = guide_legend(override.aes=list(size=1)))


