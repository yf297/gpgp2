
  double distance(double*  locs, 
                int      k, 
                int      l, 
                int     dim, 
                int     lda){
    
    double t = 0.0;

    for(int j = 0; j < dim; ++j){
        t +=  pow(locs[k + j*lda] - locs[l + j*lda],2); 
    }
    
    return( sqrt(t) );
}

  
  
  double exponential_isotropic(double*  covparms, 
                               double   d){
    
    double nugget = covparms[0] * covparms[2];
    
    if(d == 0.0){
      return(covparms[0] + nugget);
    } else {
      return(covparms[0]*exp( -(d)/covparms[1] ));
    }
  }
  
  void exponential_isotropic_mat(double*  covmat, 
                                 double*  covparms, 
                                 double*  locs, 
                                 int      bsize,
                                 int      dim,
                                 int      mb){
    
    
    for(int k = 0; k < bsize; ++k){
      for(int l = 0; l <= k; ++l){
        double d = distance(locs, k, l, dim, mb);
        covmat[k*bsize + l] = exponential_isotropic(covparms, d);
        covmat[l*bsize + k] = covmat[k*bsize + l];
      }
    }
  }


