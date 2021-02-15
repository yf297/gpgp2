source("nn_function.R")
library(GpGp)

gsize <- 100
nvec <- c(gsize,gsize)
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))
write.table(as.vector(locs), "locs.txt", row.names = F, col.names = F)

covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, nugget = 0.1)

print("simulating data...")
y <- fast_Gp_sim(covparms[c(1,2,4)], "exponential_isotropic",locs,20)
write.table(y,"y.txt", row.names = F, col.names = F)

m <- 30

print("finding nearest neighbors...")
NNarray <- find_ordered_nn(locs, m)
NNarray2 <- find_ordered_nn_brute2(locs, m)
NNarray2[is.na(NNarray2)] = 0
write.table(as.vector(t(NNarray2)),"NNarray.txt", row.names = F, col.names = F)

#print("Grouping...")
#NNlist  <- group_obs(NNarray, 2)
#mb <- get_mb(NNlist)

#print("Computing likelihoods...")
l <- system.time(c <- vecchia_meanzero_loglik(covparms, "exponential_isotropic", y, locs, NNarray))[[3]]

l
c

