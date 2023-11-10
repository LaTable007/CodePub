#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "prob.h"
#include "time.h"
#include "interface_primme.h"

int main(int argc, char *argv[])
{
  /* déclarer les variables */

  int m = 6, nev = 1;
  int n, *ia, *ja; 
  double *a;
  double *evals, *evecs;
  double tc1, tc2, tw1, tw2;
  double *datax;
  double *datay;

  /* générer le problème */
  if (prob(m, &n, &ia, &ja, &a, &datax, &datay))
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
  

  FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

  if (gnuplotPipe) {
    fprintf(gnuplotPipe, "set view 60, 30, 1, 1\n");
    fprintf(gnuplotPipe, "set ticslevel 0\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set xrange [0:2]\n");
    fprintf(gnuplotPipe, "set yrange [0:2.5]\n");
    fprintf(gnuplotPipe, "set zrange [-0.1:0.1]\n");

    for (double t = 0.0; t < 100; t += 0.1) {
        fprintf(gnuplotPipe, "splot '-' with points\n");

        // Recalcul des données pour chaque instant de temps
        for (int i = 0; i < n; i++) {
            double new_value = evecs[i] * sin(t); // Nouvelle valeur basée sur le temps
            fprintf(gnuplotPipe, "%f %f %f\n", datax[i], datay[i], new_value);
        }

        fprintf(gnuplotPipe, "e\n");
        fflush(gnuplotPipe); // Forcer l'actualisation du graphique

        usleep(100000); // Attente pour le rafraîchissement du graphique
    }

    pclose(gnuplotPipe);

  } else {
    printf("Erreur : Impossible d'ouvrir GNUplot.\n");
  }
  
  /*
  FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

  if (gnuplotPipe) {
    fprintf(gnuplotPipe, "set view 80, 30, 1, 1\n");
    fprintf(gnuplotPipe, "set ticslevel 0\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set xrange [0:2]\n");
    fprintf(gnuplotPipe, "set yrange [0:2.5]\n");
    fprintf(gnuplotPipe, "set zrange [-0.05:0.05]\n");

    fprintf(gnuplotPipe, "splot '-' with points\n");

    // Recalcul des données pour chaque instant de temps
    for (int i = 0; i < n; i++) {
      fprintf(gnuplotPipe, "%f %f %f\n", datax[i], datay[i], evecs[i]);
    }

    fprintf(gnuplotPipe, "e\n");
    pclose(gnuplotPipe);

  } 
  else {
    printf("Erreur : Impossible d'ouvrir GNUplot.\n");
  }
  */
  
  /*libérer la mémoire */
  free(ia); free(ja); free(a); free(evals); free(evecs);
  return 0;
}

