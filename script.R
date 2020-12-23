system("rm gpgp.so gpgp.o")
system("R CMD SHLIB gpgp.c")
dyn.load("gpgp.so")

library(GpGp)


vecchia_meanzero_loglik_2 <- function(covparms, y, locs, NNlist, ncores){
  
  
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
  
  ll = 0.0
  
  a <- .C("vecchia_likelihood",
          ll = as.double(ll), 
          as.double(y),
          as.double(locs), 
          as.integer(all_inds),
          as.integer(last_ind_of_block),
          as.integer(last_resp_of_block),
          as.integer(local_resp_inds),
          as.double(covparms), 
          as.integer(n),
          as.integer(nb),
          as.integer(mb),
          as.integer(dim),
          as.integer(ncores) ,
          NAOK = T)
  return(a$ll)
}


p <- 5
v <- 1:p
gsize_vec <- 60*v
m_vec <- 20*v

data1 = as.data.frame(matrix(rep(0, p*p), ncol = p))
rownames(data1) = m_vec
colnames(data1) = gsize_vec


for(i in 1:p){
  gsize <- gsize_vec[i]
  nvec <- c(gsize,gsize)
  n <- prod(nvec)
  
  covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, nugget = 0.1)
  
  x1 <- (1:nvec[1])/nvec[1]
  x2 <- (1:nvec[2])/nvec[2]
  locs <- as.matrix(expand.grid(x1,x2))
  
  y <- fast_Gp_sim(covparms[c(1,2,4)], "exponential_isotropic",locs,20)
  
  for(j in 1:p){
    m <- m_vec[j]
    NNarray <- find_ordered_nn(locs, m)
    NNlist <- group_obs(NNarray, 2)
    data1[i,j] <- system.time(vecchia_meanzero_loglik_2(covparms, y, locs, NNlist, 4))[[3]]
  }
  
}

print("time to compute mean zero exponential log-likelihood using old code")

data1= t(data1)


matplot(rownames(data1), data1, type='l', xlab='Grid Size', ylab='time (s)')
legend('topleft', inset=.05, legend=colnames(data1), 
       pch=1, horiz=TRUE,title = "Number of Neighbors",col = seq_len(ncol(data1)))
