#include <cblas-openblas.h>
#include <lapacke.h>
  
void chol(double*  a, 
          int      n){
  
    int info;
    
   LAPACKE_dpptrf(LAPACK_COL_MAJOR,'L', n, a);
    
}
  
  
void solve_l(double*  a, 
             double*  b, 
             int      n){
  
    int info;
    
    int colb[1];
    colb[0] = 1;
    
    LAPACKE_dtptrs(LAPACK_COL_MAJOR,'L', 'N', 'N', n, 1, a, b, n);

}
 
