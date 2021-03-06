#ifndef COVFUN_H
#define COVFUN_H

#include <math.h>
 
void exponential_isotropic(double*  covmat,
                           int      bsize,
                           double*  covparms, 
                           double*  locsub,
                           int      dim){
  double c0 = covparms[0];
  double c1 = covparms[1];
  double c2 = covparms[2];
	
  for(int j = 0; j < bsize; ++j){
    covmat[j*bsize + j] = c0*(1 + c2);
    for(int i = (j+1); i < bsize; ++i){
      double d = 0.0;
      for(int k = 0; k < dim; ++k){
        d += (locsub[k*bsize + j] - locsub[k*bsize + i]) * (locsub[k*bsize + j] - locsub[k*bsize + i]) ;
      }
      d = sqrt(d);
      covmat[j*bsize + i] = c0*exp(-d/c1);
    }
  }
}

#endif
