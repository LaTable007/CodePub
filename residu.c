#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void residu(int n, double *evecs, double *evals, int *ia, int *ja, double *a){
  double *res = malloc(n * sizeof(double));
  double residuNorme = 0.0;
  
  for(int i = 0; i < n; i++){
    res[i] = 0.0;
    for(int j = ia[i]; j < ia[i + 1]; j++){
      res[i] += a[j] * evecs[ja[j]];
    }
    res[i] -= evals[0] * evecs[i];
    residuNorme += (res[i] * res[i]);
  }
  residuNorme = pow(residuNorme, 1.0/2);
  printf("Norme residu = %e", residuNorme);
  free(res);
}
