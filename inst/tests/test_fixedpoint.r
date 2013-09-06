library(planar)
library(ggplot2)
library(plyr)

reflection_cell_normal <- function(lambda = NULL, k0 = 2*pi/lambda,
                                   n = list(incident=1.33, 1.5),
                                   thickness = c(30, 50)){
  
  n1 <- n[[1]]
  n2 <- n[[2]]
  k1 <- n1*k0
  k2 <- n2*k0
  l1 <- thickness[1]
  l2 <- thickness[2]
  
  r12 <- (n1 - n2) / (n1 + n2)
  r21 <- - r12
  t12 <- 2*n1 / (n1 + n2)
  t21 <- 2*n1 / (n1 + n2)
  
  r12 * exp(1i*k1*l1) + (r21*t12*t21*exp(1i*k1*l1 + 2i*k2*l2)) / 
    (1 - (r21 * exp(1i*k2*l2))^2)
}

reflection_cell <- function(lambda = NULL, k0 = 2*pi/lambda,
                            theta = NULL, q = sin(theta),
                            epsilon = list(1.33^2, 1.5^2),
                            thickness = c(30, 50)){
  
  ## define constants
  l1 <- thickness[1]
  l2 <- thickness[2]
  
  Nlambda <- length(k0)
  Nq <- length(q)
  k02 <- k0^2
  kx <- outer(k0*sqrt(epsilon[[1]] + 0i), q) # kx = q*k0
  kx2 <- kx^2
  
  kz1 <- sqrt(matrix(epsilon[[1]]*k02, nrow=Nlambda, ncol=Nq) - kx2 + 0i)  
  kz2 <- sqrt(matrix(epsilon[[2]]*k02, nrow=Nlambda, ncol=Nq) - kx2 + 0i)

  ap <- kz1 / epsilon[[1]]
  bp <- kz2 / epsilon[[2]]
  
  as <- kz1 
  bs <- kz2 
  
  r12p <-  (ap - bp) / (ap + bp)
  r12s <-  (as - bs) / (as + bs)
  r21p <-  - r12p
  r21s <-  - r12s
  t12p <-  2 * ap / (ap + bp)
  t12s <-  2 * as / (as + bs)
  t21p <-  2 * bp / (ap + bp)
  t21s <-  2 * bs / (as + bs)
  
  rp <- r12p * exp(1i*kz1*l1) + (r21p*t12p*t21p*exp(1i*kz1*l1 + 2i*kz2*l2)) / 
    (1 - (r21p * exp(1i*kz2*l2))^2)
  
  rs <- r12s * exp(1i*kz1*l1) + (r21s*t12s*t21s*exp(1i*kz1*l1 + 2i*kz2*l2)) / 
    (1 - (r21s * exp(1i*kz2*l2))^2)
  
  data.frame(q = rep(q, each=Nlambda), 
             rp = as.vector(rp), rs = as.vector(rs))
}


fixed_point <- function(lambda = NULL, k0 = 2*pi/lambda,
                        theta = NULL, q = sin(theta),
                        n = list(1.55, 1.68),
                        epsilon = lapply(n, `^`, 2),
                        thickness = c(90, 90)){
#   
#   r <- reflection_cell_normal(wavelength, 
#                     n=n, thickness=thickness)
  r <- reflection_cell(k0=k0, q=q,
                       epsilon = epsilon, 
                       thickness=thickness)
#   browser()
  rp <- r[["rp"]]
  rs <- r[["rs"]]
  leftp <-  - (rp + Conj(rp)) / (2*Mod(rp)^2)
  lefts <-  - (rs + Conj(rs)) / (2*Mod(rs)^2)
  rightp <- sqrt((rp + Conj(rp))^2 / (4*Mod(rp)^4) - 1)
  rights <- sqrt((rs + Conj(rs))^2 / (4*Mod(rs)^4) - 1)
  
  posp <- leftp + rightp
  poss <- lefts + rights
  negp <- leftp - rightp
  negs <- lefts - rights
  
  bothp <- ifelse(Mod(negp) <= 1, negp, posp)
  boths <- ifelse(Mod(negs) <= 1, negs, poss)
  
  Rp <- ifelse(Mod(bothp) <= 1, Mod(bothp)^2, NA)
  Rs <- ifelse(Mod(boths) <= 1, Mod(boths)^2, NA)
  
  data.frame(k0 = rep(k0, length(q)),
             q = r[["q"]], Rp=Rp, Rs=Rs)
}

k0lim <- 2*pi/c(100, 2000)

test <- fixed_point(k0 = seq(k0lim[1], k0lim[2], length=200), 
                    q=seq(0, 1.2, length=200),
                    thickness=c(80, 50),
                    n = c(2,3))

m <- melt(test, meas=c("Rp","Rs"))
p = ggplot(m, aes(q, k0, fill=value))+
  facet_grid(.~variable)+
  geom_tile() +
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
#   scale_fill_continuous(low="white",high="black") +
  scale_fill_continuous(low="black",high="white") +
  guides(fill="none") +
  theme(panel.background=element_rect(colour="black", size=2, fill=NA))

print(p)