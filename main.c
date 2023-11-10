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

  int m = 1, nev = 1;
  int n, *ia, *ja; 
  double *a;
  double *evals, *evecs;
  double tc1, tc2, tw1, tw2;

  /* générer le problème */
  if (prob(m, &n, &ia, &ja, &a))
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
  int i = 0;
  for(; i < n; i++){
    printf("\nVecteur propre calculée: %f\n", evecs[i]);
  }
  printf("%d\n", i);
  printf("le nombre d'élément est %d\n", n);

  double datax[] = {0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
                    0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
                    0.25, 0.5,                       1.75,
                    0.25, 0.5,                       1.75,
                    0.25, 0.5,                       1.75,
                    0.25, 0.5,                       1.75,
                    0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
                    0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
                    0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75
                    };
  double datay[] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                    0.75, 0.75, 0.75,
                    1.0, 1.0, 1.0,
                    1.25, 1.25, 1.25,
                    1.5, 1.5, 1.5,
                    1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
                    2.0, 2.0, 2.0 , 2.0, 2.0, 2.0, 2.0,
                    2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25
                    };
/*
  FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

  
  for(double t = 0.0; t < 10; t = t + 0.0001){
  
    if (gnuplotPipe) {
      fprintf(gnuplotPipe, "set view 60, 30, 1, 1\n");
      fprintf(gnuplotPipe, "set ticslevel 0\n");
      fprintf(gnuplotPipe, "set hidden3d\n");
      fprintf(gnuplotPipe, "splot '-' with points\n");

    // Écriture des données dans le fichier pour GNUplot
      for (int i = 0; i < n; i++){
       fprintf(gnuplotPipe, "%f %f %f\n", datax[i], datay[i], evecs[i] * sin(t));
      }
    
      fprintf(gnuplotPipe, "e\n");

 
      // Fermeture du fichier
      }
    else {
      printf("Erreur : Impossible d'ouvrir GNUplot.\n");
    }
  }

  pclose(gnuplotPipe);
  */
 
 

    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set view 60, 30, 1, 1\n");
        fprintf(gnuplotPipe, "set ticslevel 0\n");
        fprintf(gnuplotPipe, "set hidden3d\n");

        // Définir les plages des axes une seule fois au début
        fprintf(gnuplotPipe, "set xrange [0:2]\n");
        fprintf(gnuplotPipe, "set yrange [0:2.5]\n");
        fprintf(gnuplotPipe, "set zrange [-0.6:0.6]\n");

        for (double t = 0.0; t < 1000; t += 0.1) { // Variation de la valeur propre
            fprintf(gnuplotPipe, "splot '-' with points\n");

            // Recalculer les données pour chaque instant de temps
            for (int i = 0; i < n; i++) {
                double new_value = evecs[i] * sin(t); // Nouvelle valeur basée sur le temps
                fprintf(gnuplotPipe, "%f %f %f\n", datax[i], datay[i], new_value);
            }

            fprintf(gnuplotPipe, "e\n");
            fflush(gnuplotPipe); // Flush pour forcer l'actualisation du graphique
            if (pclose(gnuplotPipe) != -1) {
                printf("La fenêtre GNUplot a été fermée. Arrêt du programme.\n");
                break;
            }

            usleep(100000); // Pause pour laisser le temps au graphique de se rafraîchir (0.1 seconde)
        }

        pclose(gnuplotPipe);

    } else {
        printf("Erreur : Impossible d'ouvrir GNUplot.\n");
    }

  
  /*libérer la mémoire */
  free(ia); free(ja); free(a); free(evals); free(evecs);
  return 0;
}

