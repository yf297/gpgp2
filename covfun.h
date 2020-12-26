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
 
double ei(double*  covparms, 
          double   d){
   
   double nugget = covparms[0] * covparms[2];
   
   if(d == 0.0){
     return(covparms[0] + nugget);
   } else {
     return(covparms[0]*exp( -(d)/covparms[1] ));
   }
 }
 
void exponential_isotropic(double*  covmat,
                           int      bsize,
                           double*  covparms, 
                           double*  locsub,
                           int      dim){
  int ii = 0;
  for(int k = 0; k < bsize; ++k){
    for(int l = k; l < bsize; ++l){
      double d = distance(locsub, k, l, dim, bsize);
      covmat[ii + l] = ei(covparms,d);
    }
    ii += (bsize - (k+1));
  }
}


void d_exponential_isotropic(double*  dcovmat,
                             int      bsize,
                             double*  covparms, 
                             int      nparms,
                             double*  locsub,
                             int      dim){
  
  for(int k = 0; k < bsize; ++k){
    for(int l = 0; l <= bsize; ++l){
      double d = 0;
       dcovmat[0*bsize*bsize + k*bsize + l] += exp( -(d));
       dcovmat[1*bsize*bsize + k*bsize + l] += covparms[0]*exp( -(d)/covparms[1]) ;
       if(k ==l){
         dcovmat[0*bsize*bsize + k*bsize + l] += covparms[2];
         dcovmat[2*bsize*bsize + k*bsize + l] += covparms[0];
       }else{
         for(int j = 0; j < nparms; ++j){
           dcovmat[j*bsize*bsize + l*bsize + k] =dcovmat[j*bsize*bsize + k*bsize + l];
         }
       }
      }
    }
}

