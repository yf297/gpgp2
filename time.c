#include <stdio.h>
#include "src/vecchia.h"
#include <omp.h>

int main(int argc, char** argv){

    FILE *obs;
    FILE *nn;
    FILE *l;
    FILE *x;

    obs  = fopen("data/y.txt", "r");
    nn   = fopen("data/NNarray.txt", "r");
    l    = fopen("data/locs.txt", "r");
    x    = fopen("data/X.txt", "r");

    char* a = argv[1];
    char* b = argv[2];

    int n = atoi(a)*atoi(a);
    int m = atoi(b)+1;
    int p = 3;
    int dim = 2;

    double* y    = (double*) malloc(n*sizeof(double));
    double* locs = (double*) malloc(n*dim*sizeof(double));
    int* NNarray = (int*) malloc(n*m*sizeof(int));
    double* X    = (double*) malloc(n*p*sizeof(double));

    for(int  i = 0; i < n; ++i){
      fscanf(obs, "%lf", &y[i]);
    } 

    for(int i = 0; i < n*dim; ++i){
      fscanf(l, "%lf", &locs[i]); 
    }
    
    for(int i = 0; i < n*p; ++i){
      fscanf(x, "%lf", &X[i]) ;   
    }

    for(int i = 0; i < n*m; ++i){
      fscanf(nn, "%d", &NNarray[i]);
    }

    
    fclose(obs);
    fclose(nn);
    fclose(l);
    fclose(x);

    double ll = 0.0;
    double covparms[4] = {4,0.1,0.5,0.1};
    int n_cores = 4;
   
   
   double t1 = omp_get_wtime();

    vecchia_likelihood(&ll,
		       covparms, 
		       y,
               locs,
               X,
		       n,
		       dim,
		       p,
		       NNarray,
		       m,
		       n_cores);
  
   double t2 = omp_get_wtime();

   free(y);
   free(locs);
   free(NNarray);
   free(X);

   printf("likelihood:%lf\nlikelihood time:%lf\n", ll, t2-t1);
   
}
