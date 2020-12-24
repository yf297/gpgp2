#include <omp.h>
#include "util.h"
#include "linalg.h"
#include "covfun.h"
  

void vecchia_likelihood(double*  ll,
                        double* covparms,
                        double*  y, 
                        int*     n_r,
                        double*  locs,
                        int*     dim_r,
                        int*     NNarray,
                        int*     m_r, 
                        int*     ncores_r){
    
    int n = *n_r;
    int m = *m_r;
    int dim = *dim_r;
    int ncores = *ncores_r;
    
    double ySy = 0.0;
    double logdet = 0.0;
    
    double ysub[m];
    double locsub[m*dim];
    double covmat[m*m];

    for(int i = 0; i < n; ++i){
      
      int bsize = MIN(i+1, m);
      
      copy_in_ysub(ysub, bsize, y, NNarray, i*m);
      copy_in_locsub(locsub, bsize, locs, n, dim, NNarray, i*m);
    
      exponential_isotropic_mat(covmat, bsize, covparms, locsub, dim);
      
      chol(covmat, bsize);
      solve_l(covmat, ysub, bsize);
      
      ySy += pow( ysub[(bsize-1)], 2);
      logdet += 2*log(covmat[(bsize-1)*bsize]);
    }
    
    *ll = -0.5* ( n * log(2.0 * M_PI) + logdet + ySy);

}



void vecchia_grouped_likelihood(double*  ll,
                                double* covparms,
                                double*  y,
                                int*     n_r,
                                double*  locs,
                                int*     dim_r,
                                int*     all_inds,
                                int*     last_ind_of_block,
                                int*     last_resp_of_block,
                                int*     local_resp_inds,
                                int*     nb_r,
                                int*     mb_r,
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
      
      copy_in_ysub(ysub, bsize, y, all_inds, first_ind);
      copy_in_locsub(locsub, bsize, locs, n, dim, all_inds, first_ind);
      exponential_isotropic_mat(covmat, bsize, covparms, locsub, dim);
      
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

void vecchia_Linv(double* Linv,
                  double* covparms,
                  double* locs,
                  int     n,
                  int     dim,
                  int*    NNarray,
                  int     m){
    
    double locsub[m*dim];
    double covmat[m*m];

    for(int i = 0; i < n; ++i){
      
      int bsize = MIN(i+1 ,m);
      
      double one_vec[bsize];

      one_vec[bsize-1] = 0.0;

      copy_in_locsub(locsub, bsize, locs, n, dim, NNarray, i);
      exponential_isotropic_mat(covmat, bsize, covparms, locsub, dim);

      chol(covmat, bsize);
      solve_l(covmat, one_vec, bsize);

      for(int j = bsize - 1; j>=0; --j){
      
        Linv[i + (bsize-1-j)*n] = one_vec[j];
        
      }

   }             
    
}  
