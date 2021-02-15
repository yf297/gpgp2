#ifndef LINALG_H
#define LINALG_H

#include <lapacke.h>
  
void chol(double* restrict a, 
          int      n){
  
    
   LAPACKE_dpptrf(LAPACK_COL_MAJOR,'L', n, a);
    
}
  
  
void solve_l(double* restrict  a, 
             double*  b, 
             int      n){
  
        
    LAPACKE_dtptrs(LAPACK_COL_MAJOR,'L', 'N', 'N', n, 1, a, b, n);

}
 
#endif
