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
#pragma omp parallel num_threads(ncores) private(t1,t2,t3,t4,t5) 
{

  double* ysub     = (double*) malloc(m*sizeof(double) );
  double* locsub   = (double*) malloc(m*dim*sizeof(double));
  double* covmat   = (double*) malloc(mb*sizeof(double));

#pragma omp for reduction(+:ySy) reduction(+:logdet)
    for(int i = 0; i < n; ++i){
      
      int bsize = MIN(i+1, m);
      
      t1 = omp_get_wtime();
      copy_in_ungrouped(ysub,locsub, bsize, y, locs, n, dim, NNarray, i*m);
      t2 = omp_get_wtime();
      
      exponential_isotropic(covmat, bsize, covparms, locsub, dim);
      t3 = omp_get_wtime();

      chol(covmat, bsize);
      t4 = omp_get_wtime();

      solve_l(covmat, ysub, bsize);
      t5 = omp_get_wtime();

      int ii = (bsize*(bsize +1))/2;
      
       ySy += pow( ysub[(bsize-1)], 2);
       logdet += 2*log(covmat[ii-1]);
    }
    
    printf("Time to copy \n %f\n",t2-t1);
    printf("Time to fill covmat \n %f\n",t3-t2);
    printf("Time to take Cholesky \n %f\n",t4-t3);
    printf("Time to do linear solve \n %f\n",t5-t4);

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
    

    double ysub[mb];
    double locsub[dim*mb];
    double covmat[mb*mb];
    
    double t0,t00,t1,t2,t3,t4,t5;

    t0 = omp_get_wtime();
#pragma omp parallel num_threads(ncores) private(ysub,locsub,covmat,t1,t2,t3,t4,t5)
{
#pragma omp for reduction(+:ySy) reduction(+:logdet) schedule(static,ncores) 
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
      
      t1 = omp_get_wtime();
      copy_in_grouped(ysub,locsub, bsize, y, locs, n, dim, all_inds, first_ind);
      t2 = omp_get_wtime();

      exponential_isotropic(covmat, bsize, covparms, locsub, dim);
      t3 = omp_get_wtime();

      chol(covmat, bsize);
      t4 = omp_get_wtime();

      solve_l(covmat, ysub, bsize );
      t5 = omp_get_wtime();

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
    
    printf("---Times for thread number: %d\n", omp_get_thread_num());
  //printf("Time to copy in \n %f\n",t2-t1);
    printf("Time to fill covmat \n %f\n",t3-t2);
    printf("Time to take Cholesky \n %f\n",t4-t3);
  //printf("Time to do linear solve \n %f\n",t5-t4);

}    

    t00 = omp_get_wtime();
    printf("---Total time %f\n", t00-t0);
    *ll = -0.5* ( n * log(2.0 * M_PI) + logdet + ySy);
}

