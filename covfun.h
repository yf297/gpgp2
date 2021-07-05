#ifndef COVFUN_H
#define COVFUN_H

#include <math.h>
double distance(double*  locsub, 
                 int     k, 
                 int     l, 
                 int     dim, 
                 int     bsize){
    
    double t = 0.0;
    for(int j = 0; j < dim; ++j){
        t += ((locsub[k*dim + j] - locsub[l*dim + j]) * (locsub[k*dim + j] - locsub[l*dim + j])) ; 
    }

    t = sqrt(t);
    return( t );
}
 
/*double ei(double*  covparms, 
          double   d){
   
   double nugget = covparms[0] * covparms[2];
   
   if(d == 0.0){
     return(covparms[0] + nugget);
   } else {
     return(covparms[0]*exp( -(d)/covparms[1] ));
   }
 }*/
 
void exponential_isotropic(double*  covmat,
                           int      bsize,
                           double*  covparms, 
                           double*  locsub,
                           int      dim){
  double c0 = covparms[0];
  double c1 = covparms[1];
  double c2 = covparms[2];
	
  for(int k = 0; k < bsize; ++k){
    covmat[k*bsize + k] = c0*(1 + c2);
    for(int l = (k+1); l < bsize; ++l){
     double d = 0.0;
      for(int j = 0; j < dim; ++j){
        d += ((locsub[k*dim + j] - locsub[l*dim + j]) * (locsub[k*dim + j] - locsub[l*dim + j])) ; 
      }
      d = sqrt(d);
      covmat[k*bsize + l] = c0*exp(-d/c1);
    }
  }
}

#endif
