system("rm vecchia.o vecchia.so")
system("R CMD SHLIB vecchia.c")
dyn.load("vecchia.so")

library(GpGp)


vecchia_loglik_2 <- function(covparms, y, locs, NNarray, ncores){
  
  n <- length(y)
  m <- ncol(NNarray)
  dim <- ncol(locs)
  
  NNarray[is.na(NNarray)] <- 999
  inds <-c(apply(NNarray,1,sort))

  
  ll = 0.0
  ysub = rep(0,m)
  a <- .C("vecchia_likelihood",
          ll = as.double(ll), 
          as.double(covparms),
          as.double(y),
          as.integer(n),
          as.double(locs), 
          as.integer(dim),
          as.integer(inds),
          as.integer(m),
          as.integer(ncores))
  return(a$ll)
}



vecchia_grouped_loglik_2 <- function(covparms, y, locs, NNlist, ncores){


  all_inds <- as.vector(NNlist[[1]])
  last_ind_of_block <- as.vector(NNlist[[2]])
  local_resp_inds <- as.vector(NNlist[[4]])
  last_resp_of_block <- as.vector(NNlist[[5]])

  n <- length(y)
  dim <- ncol(locs)
  nb <- length(last_ind_of_block)

  x <- rep(0,nb)
  x[1] <- last_ind_of_block[1]
  for(i in 2: nb)
    x[i] <- last_ind_of_block[i] - last_ind_of_block[i-1]
  mb <- max(x)

  ll <- 0.0
 
  a<- .C("vecchia_grouped_likelihood",
         ll = as.double(ll), 
         as.double(covparms),
         as.double(y),
         as.integer(n),
         as.double(locs),
         as.integer(dim),
         as.integer(all_inds),
         as.integer(last_ind_of_block),
         as.integer(last_resp_of_block),
         as.integer(local_resp_inds),
         as.integer(nb),
         as.integer(mb),
         as.integer(ncores),
         NAOK = T)
  return(a$ll)
}




gsize <- 50
nvec <- c(gsize,gsize)
covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, nugget = 0.1)
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))
y <- fast_Gp_sim(covparms[c(1,2,4)], "exponential_isotropic",locs,20)
m = 30
NNarray <- find_ordered_nn(locs, m)
NNlist  <- group_obs(NNarray, 2)
l1 <- vecchia_grouped_loglik_2(covparms, y, locs, NNlist, 1)
l2 <- vecchia_grouped_meanzero_loglik(covparms, "exponential_isotropic", y, locs, NNlist)

l1
l2

