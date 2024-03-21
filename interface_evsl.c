#include "./EVSL_1.1.1/EVSL/include/evsl.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "./EVSL_1.1.1/EVSL/include/EVSL_config.h"
#include "./EVSL_1.1.1/EVSL/include/struct.h"

double *a;
int n, nnz, *ia, *ja, nev, ierr;

static double *vinit;
static double lmin, lmax;


int evsl(int m, int evsl_n, int* evsl_ia, int* evsl_ja, double* evsl_a, 
         int evsl_nev, double *evals, double *evecs) {
  polparams pol;
  double ecount = 1;
  double xintv[4];
  double tol = 1e-5;
  int mlan;
  ierr = EXIT_SUCCESS;

  FILE *fstats = NULL;
  n = evsl_n;
  ia = evsl_ia;
  ja = evsl_ja;
  a = evsl_a;
  nev = evsl_nev;

  int nnz = ia[n] - ia[0];
  csrMat Acsr, Bcsr;
  Acsr.ia = ia;
  Acsr.ja = ja;
  Acsr.a = a;
  Acsr.ncols = n;
  Acsr.nrows = n;
  Acsr.nnz = nnz;
  FILE *fout;

  fout = stdout;

#if EVSL_PRINT
  fstats = fout;

#else 
  fstats = fopen("output.txt","w"); 
  if (!fstats) {
    printf(" failed in opening output file in OUT/\n");
    fstats = fout;
  }

#endif 

  vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);

  EVSLStart();
  SetAMatrix(&Acsr);
  
  double t = evsl_timer();
  /*Le deuxième argument de cette fonction doit être modifié quand on 
  change la valeur de m. Pour m = 1 nous pouvons utiliser 10 et pour des 
  valeurs supérieur a 1 il est conseillé d'utiliser une grande valeur*/
  int iterr = 0;
  if(m == 1)
    iterr = 10;
  else
   iterr = 10000;
  ierr = LanTrbounds(50, iterr, tol, vinit, 1, &lmin, &lmax, fstats);//calcule las bornes maximales et minimales de l'intervalle contenant les valeurs propres

  int nev2;
  double *lam, *res;
  int *ind;
  xintv[0] = lmin;
  xintv[1] = lmax; 
  xintv[2] = lmin;
  xintv[3] = lmax;
  
  set_pol_def(&pol);
  pol.damping = 2;
  pol.thresh_int = 0.8;
  pol.thresh_ext = 0.2;
  find_pol(xintv, &pol);
  mlan = evsl_max(5 * (ecount+2), 300);
  mlan = evsl_min(mlan, n);

  ierr = ChebLanNr(xintv, mlan, tol, vinit, &pol, &nev2, &lam, &evecs, &evals, fstats);/*Calcul les valeurs propres et vecteur propres qui se trouvent dans l'intervalle fournie*/
  if (ierr) {
    printf("ChebLanTr error %d\n", ierr);
    return 1;
  }
  printf("la valeur propre minimale calculé par evsl %f\n", lam[0]);
  for(int i = 0; i < n; i++)
    printf("Composante numéro %d du vecteur propre calculé par evsl %f\n", i, evecs[i]);

  fclose(fstats);
  evsl_Free_device(vinit);

  EVSLFinish();

  return ierr; 
}
