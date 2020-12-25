setwd("~/Desktop/GpGp2")
system("rm vecchia.o vecchia.so")
system("R CMD SHLIB vecchia.c")
dyn.load("vecchia.so")
source("likelihood.R")

library(GpGp)

gsize <- 100
nvec <- c(gsize,gsize)
n  <- gsize*gsize
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, nugget = 0.1)
y <- fast_Gp_sim(covparms[c(1,2,4)], "exponential_isotropic",locs,20)

m <- 30
NNarray <- find_ordered_nn(locs, m)

NNlist  <- group_obs(NNarray, 2)
mb <- get_mb(NNlist)

l1 <- system.time(vecchia_loglik_2(covparms, y, locs, NNarray, 4))[[3]]

l2 <- system.time(vecchia_grouped_loglik_2(covparms, y, locs, NNlist, mb, 2))[[3]]

l3 <- system.time(vecchia_meanzero_loglik(covparms, "exponential_isotropic", y, locs, NNarray))[[3]]

l4 <- system.time(vecchia_grouped_meanzero_loglik(covparms, "exponential_isotropic", y, locs, NNlist))[[3]]

l1
l2
l3
l4
