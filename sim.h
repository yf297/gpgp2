#include "lmult.h"
#include "vecchia.h"
#include <imsls.h>


void fast_GP_sim_Linv(double* Linv,
                      int     n,
                      int*    NNarray,
                      int     m,
                      ){    
    
    double z[n];
 
    z = imsls_d_random_normal(n);

    L_mult(Linv, z, n, NNarray, m, y);

    }
