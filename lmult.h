#define MIN(a,b)            \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b);  \
  _a < _b ? _a : _b; })     \

  
void L_mult(double* Linv, 
            double* z,
            int     n,
            double* NNarray,
            int     m,
            double* x){

    int i;
    int j;
    int B;


    x[0] = z[0]/Linv[0 + 0*n];

    for(i=1; i<n; i++){
      
      B = MIN(i+1,m);
      x[i] = z[i];

      for(j=1; j<B; j++){
        
        x[i] -= Linv[i + j*n] * x[NNarray[i + j*n] - 1];
     
     }

      x[i] = x[i]/Linv[i];
    
    }

    return x;
}
