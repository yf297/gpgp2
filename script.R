source("likelihood.R")
library(GpGp)

gsize <- 100
nvec <- c(gsize,gsize)
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, nugget = 0.1)
y <- fast_Gp_sim(covparms[c(1,2,4)], "exponential_isotropic",locs,20)

m <- 60
NNarray <- find_ordered_nn(locs, m)

NNlist  <- group_obs(NNarray, 2)
mb <- get_mb(NNlist)

l1 <- system.time(c1 <- vecchia_loglik_2(covparms, y, locs, NNarray, 8))[[3]]

l2 <- system.time(c2 <- vecchia_grouped_loglik_2(covparms, y, locs, NNlist, mb, 8))[[3]]

l3 <- system.time(c3 <- vecchia_meanzero_loglik(covparms, "exponential_isotropic", y, locs, NNarray))[[3]]

l4 <- system.time(c4 <- vecchia_grouped_meanzero_loglik(covparms, "exponential_isotropic", y, locs, NNlist))[[3]]

l1
l2
l3
l4

