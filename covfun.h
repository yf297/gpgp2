#include <math.h>
#include <omp.h>


double distance(double*  locsub, 
                 int     k, 
                 int     l, 
                 int     dim, 
                 int     bsize){
    
    double t = 0.0;

    for(int j = 0; j < dim; ++j){
        t +=  pow(locsub[k + j*bsize] - locsub[l + j*bsize], 2); 
    }
    
    return( sqrt(t) );
}
 
 
 double exponential_isotropic(double*  covparms, 
                              double   d){
   
   double nugget = covparms[0] * covparms[2];
   
   if(d == 0.0){
     return(covparms[0] + nugget);
   } else {
     return(covparms[0]*exp( -(d)/covparms[1] ));
   }
 }


void exponential_isotropic_mat(double*  covmat,
                               int      bsize,
                               double*  covparms, 
                               double*  locsub,
                               int      dim){
  
  for(int k = 0; k < bsize; ++k){
    for(int l = k; l <= bsize; ++l){
      double d = distance(locsub, k, l, dim, bsize);
      covmat[k*bsize + l] = exponential_isotropic(covparms, d);
    }
  }
}

