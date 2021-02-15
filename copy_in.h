#ifndef COPY_IN_H
#define COPY_IN_H

void copy_in(double*  ysub,
             double*  locsub,
             int      bsize,
             double*  y,
             double*  locs,
             int      n,
             int      dim,
             int*     inds,
             int      first_ind){
  
  int ii;
  
  for(int j = 0; j < bsize; ++j){
    ii = inds[first_ind + j];
    ysub[j] = y[ii - 1];
    for(int k = 0; k < dim; ++k){
      locsub[j*dim + k]  = locs[(ii-1) + k*n];
    }
  }    
}

#endif
