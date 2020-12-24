#include "lmult.h"
#include "vecchia.h"
#include <Rmath.h>

void fast_Gp_sim_Linv(double* y,
                      int     n,
                      double* Linv,
                      int*    NNarray,
                      int     m){    

    double z[n];
    
    for(int i = 0; i < n; ++i){
        z[i] =  rnorm(0.0, 1.0);
    }
    
    L_mult(y, Linv, z, n, NNarray, m);
    
}
void fast_Gp_sim(double* y,
                 double* covparms,
                 double* locs,
                 int*    n_r,
                 int*    dim_r,
                 int*    NNarray,
                 int*    m_r){

    int n = *n_r;
    int m = *m_r;
    int dim = *dim_r;
    double Linv[n*n];
    
    vecchia_Linv(Linv, covparms, locs, n, dim, NNarray, dim);
    fast_Gp_sim_Linv(y, n, Linv, NNarray, m);
    
}


