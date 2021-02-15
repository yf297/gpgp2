#include <stdio.h>
#include "vecchia.h"

int main(){

    FILE *obs;
    FILE *nn;
    FILE *l;

    obs  = fopen("y.txt", "r");
    nn   = fopen("NNarray.txt", "r");
    l    = fopen("locs.txt", "r");

    double y[100*100];
    double locs[100*100*2];
    int NNarray[100*100*31];
    

    for(int  i = 0; i < (100*100); ++i){
      fscanf(obs, "%lf", &y[i]);
    } 

    for(int i = 0; i < (100*100*2); ++i){
      fscanf(l, "%lf", &locs[i]); 
    }
    
    for(int i = 0; i < (100*100*31); ++i){
      fscanf(nn, "%d", &NNarray[i]);
    }

    
    fclose(obs);
    fclose(nn);
    fclose(l);


    double ll = 0.0;
    double covparms[4] = {4,0.1,0.5,0.1};
    int n = 100*100;
    int dim  = 2;
    int m = 31;
    int n_cores = 4;

    double t1 = omp_get_wtime();
    vecchia_likelihood(&ll,
		       covparms, 
		       y,
		       n,
		       locs,
		       dim,
		       NNarray,
		       m,
		       n_cores);
    double t2 = omp_get_wtime();

    printf("likelihood: %lf\nlikelihood time: %lf \n", ll, t2-t1);

}
