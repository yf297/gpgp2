library(Rcpp)
#sourceCpp("~/GpGp/src/RcppExports")
sourceCpp("~/GpGp/src/matrix_mult.cpp", verbose = F)
sourceCpp("~/GpGp/src/vecchia_loglik_grad_info.cpp", verbose = F)

source("~/GpGp/R/nearest_neighbor_functions.R")
source("~/GpGp/R/ordering_functions.R")
source("~/GpGp/R/fit_model.R")
source("~/GpGp/R/simulation_functions.R")
source("nn_function.R")

gsize <- 100
nvec <- c(gsize,gsize)
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))
write.table(as.vector(locs), "data/locs.txt", row.names = F, col.names = F)

covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, nugget = 0.1)

print("simulating data...")
y <- fast_Gp_sim(covparms[c(1,2,4)], "exponential_isotropic",locs,20)
write.table(y,"data/y.txt", row.names = F, col.names = F)

m <- 80

print("finding nearest neighbors...")
NNarray <- find_ordered_nn_brute(locs, m)
NNarray2 <- find_ordered_nn_brute2(locs, m)
NNarray2[is.na(NNarray2)] = 0
write.table(as.vector(t(NNarray2)),"data/NNarray.txt", row.names = F, col.names = F)

#print("Computing likelihoods...")
t1 = Sys.time()
l <- system.time(c <- vecchia_meanzero_loglik(covparms, "exponential_isotropic", y, locs, NNarray))
t2 = Sys.time()
t2-t1
l
c
#l1 <- system.time(exponential_isotropic(covparms, locs[1:31,]))
#l1
