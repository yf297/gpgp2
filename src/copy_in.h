#ifndef COPY_IN_H
#define COPY_IN_H

void copy_in(double*  ysub,
             double*  locsub,
             double*  Xsub,
             int      bsize,
             double*  y,
             double*  locs,
             double*  X,
             int      n,
             int      dim,
             int      p,
             int*     inds,
             int      first_ind){
  
  int ii;
  
  for(int i = 0; i < bsize; ++i){
    ii = inds[first_ind + i];
    ysub[i] = y[ii - 1];
    for(int k = 0; k < dim; ++k){
      locsub[k*bsize + i]  = locs[k*n + (ii-1)];
    }
    for(int k = 0; k < p; ++k){
      Xsub[k*bsize + i] = X[k*n + (ii-1)];
    }
  }    
}

#endif
