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

  int m = 1, nev = 1;
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
    fprintf(gnuplotPipe, "set xrange [0:2]\n");
    fprintf(gnuplotPipe, "set yrange [0:2.5]\n");
    fprintf(gnuplotPipe, "set zrange [-1.0:1.0]\n");

    for (double t = 0.000; t < 100.0; t += 0.1) {
        fprintf(gnuplotPipe, "splot '-' with points\n");

        // Recalcul des données pour chaque instant de temps
        for (int i = 0; i < ne; i++) {
            if((0.75 <= datax[i]) && (datax[i] <= 1.5) && (0.75 <= datay[i]) && (datay[i] <= 1.5)){
              new_value = 0.0;
              error++;
            }
            else{
              new_value = evecs[i - error] * sin(t); // Nouvelle valeur basée sur le temps
            }
            fprintf(gnuplotPipe, "%f %f %f\n", datax[i], datay[i], new_value);
        }

        fprintf(gnuplotPipe, "e\n");
        fflush(gnuplotPipe); // Forcer l'actualisation du graphique
        usleep(100000);

    }

    pclose(gnuplotPipe);

  } else {
    printf("Erreur : Impossible d'ouvrir GNUplot.\n");
  }
  */
  

  
 /* 
   
  FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

  if (gnuplotPipe) {
    fprintf(gnuplotPipe, "set view 60, 270, 1, 1\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
      fprintf(gnuplotPipe, "set dgrid3d 50,50 qnorm 5\n");
    fprintf(gnuplotPipe, "set xrange [0:2]\n");
    fprintf(gnuplotPipe, "set yrange [0:2.5]\n");
    fprintf(gnuplotPipe, "set zrange [-1.0:1.0]\n");

    fprintf(gnuplotPipe, "splot '-' with lines\n");
  
    // Recalcul des données pour chaque instant de temps
    for (int i = 0; i < ne; i++) {
      if((0.75 <= datax[i]) && (datax[i] <= 1.5) && (0.75 <= datay[i]) && (datay[i] <= 1.5)){
        new_value = 0.0;
        error++;
      }
      else {
        new_value = evecs[i - error];
      }
      printf("%f %f %f\n", datax[i], datay[i], new_value);
      fprintf(gnuplotPipe, "%f %f %f\n", datax[i], datay[i], new_value);
    }

    fprintf(gnuplotPipe, "e\n");
    pclose(gnuplotPipe);

  } 
  else {
    printf("Erreur : Impossible d'ouvrir GNUplot.\n");
  }
  */
  

  //Méthode d'Euler progressive

  double T0 = 20;
  double D = 9.7 / (10 * 10 * 10 * 10 * 10);
  double h = 1.0;




  //Calculons f(t, u(t))
  
  double *U = malloc(n * sizeof(double));
  double *result = malloc(n * sizeof(double));

  

  for(int i = 0; i < n; i++){
    U[i] = T0;
    printf("U[%d] = %f\n", i, U[i]);
  }
  
  /*
  FILE *fp = NULL; // Ouvrir le fichier pour écrire les données
  fp = popen("gnuplot -persistent", "w");
  

  fprintf(fp, "viewmap");
  fprintf(fp, "dgrid3d");
  
  for (double k = 0.0; k < 100; k = k + 0.1){
    fprintf(fp, "splot '-' with pm3d\n");

    
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            result[i] += (-D) * U[j] * a[j]; 
        }
        fprintf(fp, "%f %f %f\n", datax[i], datay[i], U[i]);
        U[i] = U[i] + (h * result[i]);
    }
    fprintf(fp, "e\n");
    fflush(fp); 
  }
  fclose(fp);
  */
  
  FILE *fp = NULL; // Ouvrir le fichier pour écrire les données
  fp = popen("gnuplot -persist", "w");
  if(fp){
    //fprintf(fp, "unset key\n");
    //fprintf(fp, "set view map\n");
    //fprintf(fp, "set pm3d interpolate 10,10\n");
    //fprintf(fp, "set cbrange [0:25]\n");
    
    
    for (double k = 0; k < 2; k += 1) {
      //fprintf(fp, "splot '-' with pm3d\n");
      printf("Itération numéro %f\n\n", k);
     
      for (int i = 0; i < ne; i++) {
        result[i] = 0.0;
        if((0.75 <= datax[i]) && (datax[i] <= 1.5) && (0.75 <= datay[i]) && (datay[i] <= 1.5)){
          error++;
          //printf("%f %f %f %d\n", datay[i], datax[i], 0.0, i);
          //fprintf(fp, "%f %f %f\n", datay[i], datax[i], 0.0);
          continue;
        }
        matvec(U, result, &nev, n, ia, ja, a);
         
        for(int j = 0; j < n; j++){
          result[j - error] = -D * result[j - error];
          U[i - error] = U[i - error] + (h * result[i - error]);
          printf("%f / %f\n", result[j - error], U[i - error]);
        }

        //printf("%f %f %f %d\n", datay[i], datax[i], U[i], i);
        //fprintf(fp, "%f %f %f\n", datay[i], datax[i], U[i]);
          
        if((i + 1) % nx == 0){ 
          //printf("\n");
          //fprintf(fp, "\n");
        }
      }
      error = 0.0;
      //printf("e\n");
      //fprintf(fp, "e\n");
      //printf("\n");




      // Marquer la fin des données pour chaque ensemble
      
      //fflush(fp); // Forcer l'actualisation du graphique
      //usleep(100000);
    }
    pclose(fp);
  }
  
  


  /*libérer la mémoire */
  free(ia); free(ja); free(a); free(evals); free(evecs); free(datax); free(datay); free(U); free(result);
  return 0;
}

