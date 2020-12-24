#include "util.h"
  
void L_mult(double* y,
            double* Linv, 
            double* z,
            int     n,
            int* NNarray,
            int     m){

    int i;
    int j;
    int B;

    y[0] = z[0]/Linv[0 + 0*n];

    for(i=1; i<n; i++){
      
      B = MIN(i+1,m);
      y[i] = z[i];

      for(j=1; j<B; j++){
        
        y[i] += Linv[i + j*n] * y[NNarray[i + j*n] - 1];
     
     }

      y[i] = y[i]/Linv[i];
    
    }
}
