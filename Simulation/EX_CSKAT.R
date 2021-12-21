.libPaths("~/tmp")
# install.packages("MiRKAT")

install.packages("~/MiRKAT_1.1.1.tar.gz", repos=NULL, type="source")

library(MiRKAT)
Y <- rnorm(100)
Z <- matrix(rnorm(200), 100, 2)
ID <- gl(20, 5)
G <- matrix(rbinom(1000, 2, 0.05), 100, 10)
K <- G %*% t(G)
CSKAT(formula.H0 = Y ~ Z + (1 | ID), Ks = K)


