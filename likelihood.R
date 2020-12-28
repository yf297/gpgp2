dyn.load("vecchia")

vecchia_loglik_2 <- function(covparms, y, locs, NNarray, ncores){
  
  n <- length(y)
  m <- ncol(NNarray)
  dim <- ncol(locs)
  ll <- 0.0
  nparms <- length(covparms)
  
  a <- .C("vecchia_likelihood",
          ll = as.double(ll),
          as.double(covparms),
          as.integer(nparms),
          as.double(y),
          as.integer(n),
          as.double(locs),
          as.integer(dim),
          as.integer(t(NNarray)),
          as.integer(m),
          as.integer(ncores),
          NAOK = TRUE)
  return(a$ll)
}



vecchia_grouped_loglik_2 <- function(covparms, y, locs, NNlist, mb, ncores){
  
  all_inds <- as.vector(NNlist[[1]])
  last_ind_of_block <- as.vector(NNlist[[2]])
  local_resp_inds <- as.vector(NNlist[[4]])
  last_resp_of_block <- as.vector(NNlist[[5]])
  
  n <- length(y)
  dim <- ncol(locs)
  nb <- length(last_ind_of_block)
  
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


get_mb <- function(NNlist){
  
  all_inds <- as.vector(NNlist[[1]])
  last_ind_of_block <- as.vector(NNlist[[2]])
  local_resp_inds <- as.vector(NNlist[[4]])
  last_resp_of_block <- as.vector(NNlist[[5]])
  nb <- length(last_ind_of_block)
  x <- rep(0,nb)
  x[1] <- last_ind_of_block[1]
  
  for(i in 2: nb)
    x[i] <- last_ind_of_block[i] - last_ind_of_block[i-1]
  mb <- max(x)
  
  return(mb)
}



