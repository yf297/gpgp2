#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>


  
void chol(double*  a, 
          int      n){
  
    int info;
    
    F77_CALL(dpptrf)("L", &n, a, &info);
    
}
  
  
void solve_l(double*  a, 
             double*  b, 
             int      n){
  
    int info;
    
    int colb[1];
    colb[0] = 1;
    
    F77_CALL(dtptrs)("L", "N", "N", &n, colb, a, b, &n, &info);
}
 