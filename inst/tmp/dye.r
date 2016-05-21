
# bill paper

epsr <- read.table("epsr.txt", head=TRUE)
epsi <- read.table("epsi.txt", head=TRUE)

wvl <- seq(300, 800)
bill <- data.frame(wavelength = wvl, 
                                epsilon = predict(smooth.spline(epsr$x,epsr$y), wvl)$y + 
                                  1i*predict(smooth.spline(epsi$x,epsi$y), wvl)$y)

saveRDS(dye, file="bill.rds")
# write.table(dye, "dye.txt", row.names = FALSE)
dye <- readRDS("bill.rds")

dye_fun <- function(wavelength){
  eps_r <- predict(smooth.spline(dye$wavelength, Re(dye$epsilon)), wavelength)$y
  eps_i <- predict(smooth.spline(dye$wavelength, Im(dye$epsilon)), wavelength)$y
  data.frame(wavelength = wavelength, epsilon = eps_r + 1i*eps_i, real=eps_r, imag=eps_i)
}

test <- dye_fun(seq(300, 800))


m <- reshape2::melt(test, id="wavelength", meas=c("real", "imag"))
ggplot(m, aes(wavelength, value))+ facet_wrap(~variable, scales="free", ncol=1) + geom_line()




# cc paper

n <- read.table("n.txt", head=TRUE)
k <- read.table("k.txt", head=TRUE)

wvl <- seq(300, 800)
dye <- dplyr::mutate(data.frame(wavelength = wvl, 
                                n = predict(smooth.spline(n$x,n$y), wvl)$y + 1i*predict(smooth.spline(k$x,k$y), wvl)$y), 
                     epsilon = n^2)

saveRDS(dye, file="dye.rds")
# write.table(dye, "dye.txt", row.names = FALSE)
dye <- readRDS("dye.rds")

dye_fun <- function(wavelength){
  eps_r <- predict(smooth.spline(dye$wavelength, Re(dye$epsilon)), wavelength)$y
  eps_i <- predict(smooth.spline(dye$wavelength, Im(dye$epsilon)), wavelength)$y
  data.frame(wavelength = wavelength, epsilon = eps_r + 1i*eps_i, real=eps_r, imag=eps_i)
}

test <- dye_fun(seq(300, 800))


m <- reshape2::melt(test, id="wavelength", meas=c("real", "imag"))
ggplot(m, aes(wavelength, value))+ facet_wrap(~variable, scales="free", ncol=1) + geom_line()

