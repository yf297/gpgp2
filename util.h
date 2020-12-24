#define MIN(a,b)            \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b);  \
  _a < _b ? _a : _b; })     \


void copy_in_ysub(double*   ysub,
                  int       bsize,
                  double*   y,
                  int*      inds,
                  int       first_ind){
    
    for(int j = 0; j < bsize; ++j){
      ysub[j] = y[ inds[first_ind + j] - 1];
    }  
}
 
void copy_in_locsub(double* locsub,
                    int     bsize,
                    double* locs,
                    int     n,
                    double  dim,
                    int*    inds,
                    int     first_ind){
    
    for(int j = 0; j < bsize; ++j){
      for(int k = 0; k < dim; ++k){
        locsub[j + k*bsize] = locs[ (inds[first_ind + j] - 1) + k*n];
      }
    }      
}
                    
