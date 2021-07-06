#ifndef LINALG_H
#define LINALG_H

#include <cblas.h>
#include <lapacke.h>
  
void chol(double*  a, 
          int      n){
  
    
   LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L', n, a, n);
    
}
  
  
void solve(double* a, 
             double*  b, 
             int      n){
  
        
    LAPACKE_dtrtrs(LAPACK_COL_MAJOR,'L', 'N', 'N', n, 1, a,n, b, n);

}

void solve_t(double* a, 
             double*  b, 
             int      n){
  
        
    LAPACKE_dtrtrs(LAPACK_COL_MAJOR,'L', 'T', 'N', n, 1, a,n, b, n);

}

double dot(double* a,
         double* b,
         int     n){
    
   return cblas_ddot(n, a, 1, b, 1);
    
    }

#endif
