library(planar)

lambda <- 500
epsilon <- list(1.5^2, -20+1i, 1.33^2)
thickness <- c(0,50,0)

              
Mtot1 <- dipole(lambda=lambda, epsilon=epsilon, thickness=thickness,
                Nquadrature1 = 0, Nquadrature2 = 0,
                Nquadrature3 = 0, qcut = NULL, GL=FALSE)

Mtot2 <- dipole(lambda=lambda, epsilon=epsilon, thickness=thickness,
                Nquadrature1 = 200, Nquadrature2 = 500,
                Nquadrature3 = 200, qcut = NULL, GL=TRUE)

Mtot3 <- dipole_direct(lambda=lambda, epsilon=epsilon, thickness=thickness,
                Nquadrature1 = 500, Nquadrature2 = 1e3,
                Nquadrature3 = 500, qcut = NULL)


