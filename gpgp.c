#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include <omp.h>
#include <R_ext/Lapack.h>



#define MIN(a,b)            \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b);  \
  _a < _b ? _a : _b; })     \
    
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


  
  
void chol(double*  a, 
          int      n){
  
    int info;

    
    F77_CALL(dpotrf)("L", &n, a, &n, &info);
    
    for (int j = 0; j < n; ++j) {
      for (int i = 1; i < n; ++i) {
        a[j*n + (i-j)] = a[j*n + i];
      }
    }
}
  
  
void solve_l(double*  a, 
             double*  b, 
             int      n){
  
    int info;
    
    int sub[1];
    sub[0] = n-1;
    
    int colb[1];
    colb[0] = 1;
    
    F77_CALL(dtbtrs)("L", "N", "N", &n, sub, colb, a, &n, b, &n, &info);
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



  
void copy_in(  double*   ysub, 
               double*   locsub,
               double*   y,
               double*   locs,
               int*      all_inds,
               int       first_ind,
               int       bsize,
               int       n,
               int       mb,
               int       dim)
{
    

      for(int j = 0; j < bsize; ++j){
        ysub[j] = y[ all_inds[first_ind + j] - 1  ];
        for(int k=0; k< dim; ++k){
          locsub[j + k*mb] = locs[all_inds[first_ind + j] - 1 + k*n];
        }
      }
      
}
  
  

void vecchia_likelihood(  double*  ll, 
                          double*  y, 
                          double*  locs,
                          int*     all_inds,
                          int*     last_ind_of_block,
                          int*     last_resp_of_block,
                          int*     local_resp_inds,
                          double*  covparms,
                          int*     n_r,
                          int*     nb_r,
                          int*     mb_r,
                          int*     dim_r,
                          int*     ncores_r){
    int n = *n_r;
    int nb = *nb_r;
    int mb = *mb_r;
    int dim = *dim_r;
    int ncores = *ncores_r;
    
    double ySy = 0.0;
    double logdet = 0.0;
    

    double ysub[mb];
    double locsub[dim*mb];
    double covmat[mb*mb];

#pragma omp parallel for num_threads(ncores) private(covmat,ysub,locsub) reduction(+:ySy) reduction(+:logdet) schedule(static)
    for(int i = 0; i < nb; ++i){
      
      int first_ind;
      int last_ind;
      int bsize;
      
      if(i==0){ first_ind = 0; } else {first_ind = last_ind_of_block[i-1] + 1 - 1; }
      last_ind = last_ind_of_block[i] - 1;
      bsize= last_ind - first_ind + 1;
      
      int first_resp;
      int last_resp;
      int rsize;
      
      if(i == 0){ first_resp = 0; } else {first_resp = last_resp_of_block[i-1]+1-1; }
      last_resp = last_resp_of_block[i]-1;
      rsize = last_resp - first_resp + 1;
      
      copy_in(ysub,locsub, y, locs, all_inds, first_ind,bsize,n, mb,dim);
      
      exponential_isotropic_mat(covmat,covparms,locsub,bsize, dim,mb);
      
      chol(covmat, bsize);
      solve_l(covmat, ysub, bsize );
      
      int wr;
      for(int j = 0; j < rsize; ++j){
        wr  =  local_resp_inds[first_resp+j] - 1;
        logdet += 2.0*log( covmat[ wr * bsize ] ) ; 
        ySy +=  pow( ysub[wr], 2);
      }
      
    }
    
    *ll = -0.5* ( n * log(2.0 * M_PI) + logdet + ySy);
    
  }
