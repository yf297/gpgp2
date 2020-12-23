#include <omp.h>
#include "linalg.h"
#include "covfun.h"
  
 void copy_in( double*   ysub, 
               double*   locsub,
               double*   y,
               double*   locs,
               int*      all_inds,
               int       first_ind,
               int       bsize,
               int       n,
               int       mb,
               int       dim){
    
    for(int j = 0; j < bsize; ++j){
      ysub[j] = y[ all_inds[first_ind + j] - 1  ];
      for(int k=0; k< dim; ++k){
        locsub[j + k*mb] = locs[all_inds[first_ind + j] - 1 + k*n];
      }
    }  
 }
  
  

 void vecchia_likelihood( double*  ll, 
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
