#ifndef LINALG_H
#define LINALG_H

#include <lapacke.h>
  
void chol(double*  a, 
          int      n){
  
    
   LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L', n, a, n);
    
}
  
  
void solve_l(double* a, 
             double*  b, 
             int      n){
  
        
    LAPACKE_dtrtrs(LAPACK_COL_MAJOR,'L', 'N', 'N', n, 1, a,n, b, n);

}
 
#endif
