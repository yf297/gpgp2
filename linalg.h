
#include <R_ext/Lapack.h>


  
void chol(double*  a, 
          int      n){
  
    int info;
    
    F77_CALL(dpotrf)("L", &n, a, &n, &info);
    
    for (int j = 0; j < n; ++j) {
      for (int i = 1; i < n; ++i) {
        a[j*n + (i-j)] = a[j*n + i];
      }
    }
}
  
  
void solve_l(double*  a, 
             double*  b, 
             int      n){
  
    int info;
    
    int sub[1];
    sub[0] = n-1;
    
    int colb[1];
    colb[0] = 1;
    
    F77_CALL(dtbtrs)("L", "N", "N", &n, sub, colb, a, &n, b, &n, &info);
}
 
