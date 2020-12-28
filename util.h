#define MIN(a,b)            \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b);  \
  _a < _b ? _a : _b; })     \
  

void copy_in_ungrouped(double* restrict  ysub,
                       double* restrict  locsub,
                       int       bsize,
                       double* restrict   y,
                       double* restrict  locs,
                       int       n,
                       int       dim,
                       int*      inds,
                       int       first_ind){
  
  int ii;
  
  for(int j = bsize-1; j >= 0; --j){
    ii = inds[first_ind + j];
    ysub[bsize-1-j] = y[ii - 1];
  }
  
  for(int k = 0 ; k < dim; ++k){
    for(int j = bsize-1; j >= 0; j--){
      ii = inds[first_ind + j];
      locsub[bsize-1-j + k*bsize] = locs[(ii - 1) + k*n];
    }  
  }
  
}

void copy_in_grouped(double* restrict   ysub,
                     double* restrict  locsub,
                     int       bsize,
                     double* restrict  y,
                     double* restrict  locs,
                     int       n,
                     int       dim,
                     int*      inds,
                     int       first_ind){
  
  int ii;
  
  for(int j = 0; j < bsize; ++j){
    ii = inds[first_ind + j];
    ysub[j] = y[ii - 1];
  }
  
  for(int k = 0 ; k < dim; ++k){
    for(int j = 0; j < bsize; ++j){
      ii = inds[first_ind + j];
      locsub[j + k*bsize] = locs[(ii - 1) + k*n];
    }  
  }
  
}



void copy_in_locsub(double* locsub,
                    int     bsize,
                    double* locs,
                    int     n,
                    double  dim,
                    int*    inds,
                    int     first_ind){
  
  for(int k = 0; k < dim; ++k){
    for(int j = 0; j < bsize; ++j){
      locsub[j + k*bsize] = locs[ (inds[first_ind + j] - 1) + k*n];
    }
  }      
}        
