#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "linalg.h"
#include "covfun.h"
#include "util.h"

void vecchia_likelihood(double*  ll,
                        double*  covparms,
                        int*     nparms_r,
                        double*  y,
                        int*     n_r,
                        double*  locs,
                        int*     dim_r,
                        int*     NNarray,
                        int*     m_r,
                        int*     ncores_r){


    double t1,t2,t3,t4,t5;


    int n = *n_r;
    int m = *m_r;
    int dim = *dim_r;
    int ncores = *ncores_r;
    
    double ySy = 0.0;
    double logdet = 0.0;
    
    int mb = m*(m+1)/2;
#pragma omp parallel num_threads(ncores) 
{

  double* ysub     = (double*) malloc(m*sizeof(double) );
  double* locsub   = (double*) malloc(m*dim*sizeof(double));
  double* covmat   = (double*) malloc(mb*sizeof(double));

#pragma omp for reduction(+:ySy) reduction(+:logdet) schedule(static)
    for(int i = 0; i < n; ++i){
      
      int bsize = MIN(i+1, m);
      
      copy_in_ungrouped(ysub,locsub, bsize, y, locs, n, dim, NNarray, i*m);
      
      exponential_isotropic(covmat, bsize, covparms, locsub, dim);

      chol(covmat, bsize);

      solve_l(covmat, ysub, bsize);

      int ii = (bsize*(bsize +1))/2;
      
      ySy += pow( ysub[(bsize-1)], 2);
      logdet += 2*log(covmat[ii-1]);
    }
    

}
    *ll = -0.5* (n * log(2.0 * M_PI) + logdet + ySy);
     

}



void vecchia_grouped_likelihood(double*  ll,
                                double*  covparms,
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
    


#pragma omp parallel num_threads(ncores)
{

      double* ysub    = (double*) malloc(mb*sizeof(double));
      double* locsub  = (double*) malloc(mb*mb*sizeof(double));
      double* covmat  = (double*) malloc(mb*mb*sizeof(double));


#pragma omp for reduction(+:ySy) reduction(+:logdet) schedule(static) 
    for(int i = 0; i < nb; ++i){
      
      int first_ind;
      int last_ind;
      int bsize;
      
      if(i==0){ first_ind = 0; } else {first_ind = last_ind_of_block[i-1] + 1 - 1; }
      last_ind = last_ind_of_block[i] - 1;
      bsize = last_ind - first_ind + 1;
      
      int first_resp;
      int last_resp;
      int rsize;
      
      if(i == 0){ first_resp = 0; } else {first_resp = last_resp_of_block[i-1]+1-1; }
      last_resp = last_resp_of_block[i]-1;
      rsize = last_resp - first_resp + 1;


      copy_in_grouped(ysub,locsub, bsize, y, locs, n, dim, all_inds, first_ind);  
     
      exponential_isotropic(covmat, bsize, covparms, locsub, dim);
      
      chol(covmat, bsize);

      solve_l(covmat, ysub, bsize );
      
      int wr;
      int ii;
      
      for(int j = 0; j < rsize; ++j){
      ii = 0;
        wr  =  local_resp_inds[first_resp + j] - 1;
        for(int s = 0; s < wr; ++s){
          ii += (bsize - s);
        }
        logdet += 2.0*log( covmat[ ii ] ) ; 
        ySy +=  pow( ysub[wr], 2);
      }
    }
    

}    

    *ll = -0.5* ( n * log(2.0 * M_PI) + logdet + ySy);
}

