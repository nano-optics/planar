## SCARY normalisation factors when impedance mismatch between incident and outgoing media
## this seems to work now

library(planar)
stackp <- list(epsilon=list(1.2, rep(1,2), 1.5),
               wavelength=c(500, 550), 
               thickness=c(0, 200, 0),
               angle=0*pi/180, polarisation='p')
# 
# stackp <- list(epsilon=list(1.5^2, (0.2+1i)^2, 1.3^2),
#                wavelength=c(500), 
#                thickness=c(0, 50, 0),
#                angle=0*pi/180, polarisation='p')

stacks <- modifyList(stackp, list(polarisation="s"))

simul <- function(stack){
  
  res <- data.frame(do.call(multilayer, stack))
  #   res <- data.frame(do.call(multilayercpp, stack))
    #res <- data.frame(do.call(recursive_fresnelcpp, stack))
  res$T2 <- 1 - res$R
  res
}

test <- simul(stackp)
test2 <- simul(stacks)

with(test, matplot(k0, cbind(T,T2), t="p",col=1))
# with(test, points(k0, T2, t="p"))
with(test2, points(k0, T, t="p",col="red", pch="+"))
