#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "prob.h"
#include "time.h"
#include "interface_primme.h"
#include "matvec.h"

int main(int argc, char *argv[])
{
  /* déclarer les variables */

  int m = 2, nev = 1;
  int n, *ia, *ja; 
  double *a;
  double *evals, *evecs;
  double tc1, tc2, tw1, tw2;
  double *datax;
  double *datay;
  int ne;
  int nx;

  /* générer le problème */
  if (prob(m, &n, &ia, &ja, &a, &datax, &datay, &ne, &nx))
     return 1;

  printf("\nPROBLÈME: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
  
  /* allouer la mémoire pour vecteurs & valeurs propres */
  evals = malloc(nev * sizeof(double));
  evecs = malloc(nev * n * sizeof(double));

  if (evals == NULL || evecs == NULL) {
     printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs propres\n\n");
      return 1;
  }

  /* primme - résolution */
  tc1 = mytimer_cpu(); tw1 = mytimer_wall();
  if(primme(n, ia, ja, a, nev, evals, evecs))
     return 1;
  tc2 = mytimer_cpu(); tw2 = mytimer_wall();

  /* temps de solution */
  printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1);
  printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1);
  printf("\nValeur propre minimale calculée: %5.1f\n",evals[0]);
  double new_value;
  int error = 0.0;
  
 /* 
  FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

  if (gnuplotPipe) {
    fprintf(gnuplotPipe, "set view 60, 30, 1, 1\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set pm3d\n");
    //fprintf(gnuplotPipe, "set cbrange [-0.35:0.35]\n");
    //fprintf(gnuplotPipe, "set zrange [-0.35:0.35]\n");



      
    // Recalcul des données pour chaque instant de temps
    for(double t = 0.0; t < 0.1; t = t + 0.1){
      fprintf(gnuplotPipe, "splot '-' with pm3d\n");

    
      for (int i = 0; i < ne; i++) {
        if((3.0 <= datax[i]) && (datax[i] <= 6.0) && (3.0 <= datay[i]) && (datay[i] <= 6.0) ){
          //printf("%f %f %f\n", datay[i], datax[i], 0.0);
          fprintf(gnuplotPipe, "%f %f %f\n", datay[i], datax[i], 0.0);
          error++;
        }
        else if((datay[i] == 0.0) || (datay[i] == 10.0) || (datax[i] == 0.0) || (datax[i] == 8.0)){
          //printf("%f %f %f\n", datay[i], datax[i], 0.0);
          fprintf(gnuplotPipe, "%f %f %f\n", datay[i], datax[i], 0.0);
          error++;
        }

        else {
          //printf("%f %f %d %f\n", datay[i], datax[i], error, evecs[i - error] * sinf(t * evecs[0]));
          //fprintf(gnuplotPipe, "%f %f %f\n", datay[i], datax[i], evecs[i - error] * sinf(t * evals[0]) * m);
          fprintf(gnuplotPipe, "%f %f %f\n", datay[i], datax[i], evecs[i - error] * m);
        }
        
        if((i + 1) % (nx + 2) == 0){ 
          //printf("\n");
          fprintf(gnuplotPipe, "\n");
        }
      }
      //printf("e\n");
      fprintf(gnuplotPipe, "e\n");
  
      
      fflush(gnuplotPipe);
      usleep(100000);
      error = 0.0;
    }
    pclose(gnuplotPipe);
  }

  else {
    printf("Erreur : Impossible d'ouvrir GNUplot.\n");
  }
  */
  

  //Méthode d'Euler progressive

  double T0 = 20;
  double D = 9.7 / (10 * 10 * 10 * 10 * 10);
  double h = 50.0;




  //Calculons f(t, u(t))
  
  double *U = malloc(n * sizeof(double));
  double *result = malloc(n * sizeof(double));

  

  for(int i = 0; i < n; i++){
    U[i] = T0;
    //printf("U[%d] = %f\n", i, U[i]);
  }
 
  FILE *fp = NULL; // Ouvrir le fichier pour écrire les données
  fp = popen("gnuplot -persist", "w");
  if(fp){
    fprintf(fp, "unset key\n");
    fprintf(fp, "set view map\n");
    fprintf(fp, "set pm3d interpolate 10,10\n");
    fprintf(fp, "set cbrange [0:20]\n");
    
    
    for (double k = 0; k < 2000; k++) {
      fprintf(fp, "splot '-' with pm3d\n");
      //printf("Itération numéro %f\n\n", k);


      for(int i = 0; i < ne; i++){
        if((3.0 <= datax[i]) && (datax[i] <= 6.0) && (3.0 <= datay[i]) && (datay[i] <= 6.0)){
          error++;
          //printf("%f %f %f %d\n", datay[i], datax[i], 0.0, i);
          fprintf(fp, "%f %f %f\n", datay[i], datax[i], 0.0);
          continue;
        }
        else if((datay[i] == 0.0) || (datay[i] == 10.0) || (datax[i] == 0.0) || (datax[i] == 8.0)){
          //printf("%f %f %f %d\n", datay[i], datax[i], 0.0, i);
          fprintf(fp, "%f %f %f\n", datay[i], datax[i], 0.0);
          error++;
        }
        else{
          fprintf(fp, "%f %f %f\n", datay[i], datax[i], U[i - error]);
          //printf("%f %f %f %d\n", datay[i], datax[i], U[i - error], i);
        }
          
        if((i + 1) % (nx + 2) == 0){ 
          //printf("\n");
          fprintf(fp, "\n");
        }
      }

  
      matvec(U, result, &nev, n, ia, ja, a);
         
      for(int j = 0; j < n; j++){
        result[j] = -D * result[j];
        U[j] = U[j] + (h * result[j]);
        //printf("%f / %f\n", result[j - error], U[j - error]);
      }
      error = 0.0;

      //printf("e\n");
      fprintf(fp, "e\n");
      //printf("\n");

      
      fflush(fp); // Forcer l'actualisation du graphique
      usleep(100000);
    }

    pclose(fp);
  }
  

 


  /*libérer la mémoire */
  free(ia); free(ja); free(a); free(evals); free(evecs); free(datax); free(datay); free(U); free(result);
  return 0;
}

