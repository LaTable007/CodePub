#include "interface_evsl.h"
#include "interface_primme.h"
#include "matvec.h"
#include "primme/PRIMMESRC/COMMONSRC/primme.h"
#include "prob.h"
#include "residu.h"
#include "time.h"
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
  /* déclarer les variables */
  int m = 30, nev = 1;
  int n, *ia, *ja;
  double *a;
  double *evals, *evecs;
  double tc1, tc2, tc3, tc4, tw1, tw2, tw3, tw4;
  double *datax;
  double *datay;
  int ne;
  int nx;

  /* initialiser le générateur de nombres aléatoires */
  /* générer le problème */
  if (prob(m, &n, &ia, &ja, &a, &datax, &datay, &ne, &nx))
    return 1;

  printf("\nPROBLÈME: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n]);

  /* allouer la mémoire pour vecteurs & valeurs propres */
  evals = malloc(nev * sizeof(double));
  evecs = malloc(nev * n * sizeof(double));

  if (evals == NULL || evecs == NULL) {
    printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs "
           "propres\n\n");
    return 1;
  }

  /* primme - résolution */
  tc1 = mytimer_cpu();
  tw1 = mytimer_wall();
  if (primme(n, ia, ja, a, nev, evals, evecs, primme_smallest))
    return 1;
  tc2 = mytimer_cpu();
  tw2 = mytimer_wall();

  /* temps de solution */
  printf("\nTemps de solution (CPU): %5.1f sec", tc2 - tc1);
  printf("\nTemps de solution (horloge): %5.1f sec \n", tw2 - tw1);
  printf("\nValeur propre minimale calculée: %5.1f\n\n", evals[0]);

  /*PARTIE RESIDU*/

  tc3 = mytimer_cpu();
  tw3 = mytimer_wall();
  residu(n, evecs, evals, ia, ja, a);
  tc4 = mytimer_cpu();
  tw4 = mytimer_wall();

  printf("\nTemps de solution (CPU): %e sec", tc4 - tc3);
  printf("\nTemps de solution (horloge): %e sec \n\n", tw4 - tw3);

  int error = 0.0;

  /*Affichage vecteur propre*/

  FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

  if (gnuplotPipe) {
    fprintf(gnuplotPipe, "set view 60, 210, 1, 1\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set pm3d\n");

    // Recalcul des données pour chaque instant de temps
    for (int t = 0; t < 1; t++) {
      fprintf(gnuplotPipe, "splot '-' with pm3d\n");

      for (int i = 0; i < ne; i++) {
        if ((3.0 <= datax[i]) && (datax[i] <= 6.0) && (3.0 <= datay[i]) &&
            (datay[i] <= 6.0)) { /*On vérifie si on se trouve dans le trou ou
                                    sur le bord du trou*/
          fprintf(gnuplotPipe, "%f %f %f\n", datay[i], datax[i], 0.0);
          error++;
        } else if ((datay[i] == 0.0) || (datay[i] == 10.0) ||
                   (datax[i] == 0.0) ||
                   (datax[i] == 8.0)) { /*On vérifie si on se trouve sur le bord
                                           de la grile*/
          fprintf(gnuplotPipe, "%f %f %f\n", datay[i], datax[i], 0.0);
          error++;
        }

        else {
          fprintf(gnuplotPipe, "%f %f %f\n", datay[i], datax[i],
                  evecs[i - error + (n * (nev - 1))] * m);
        }

        if ((i + 1) % (nx + 2) == 0) {
          fprintf(gnuplotPipe, "\n");
        }
      }
      fprintf(gnuplotPipe, "e\n");

      fflush(gnuplotPipe);
      error = 0.0;
    }
    pclose(gnuplotPipe);
  }

  else {
    printf("Erreur : Impossible d'ouvrir GNUplot.\n");
  }

  // Méthode d'Euler progressive

  double T0 = 20;
  double D = 9.7 / (10 * 10 * 10 * 10 * 10);

  double *U = malloc(n * sizeof(double));
  double *result = malloc(n * sizeof(double));

  double *pvals = malloc(nev * sizeof(double));
  double *pvecs = malloc(nev * n * sizeof(double));
  double *Da = malloc(ia[n] * sizeof(double));

  for (int i = 0; i < ia[n]; i++) {
    Da[i] = D * a[i];
  }

  if (primme(n, ia, ja, Da, nev, pvals, pvecs, primme_largest))
    return 1;

  printf("Valeur propre maximale calculée: %f\n", pvals[0]);

  printf("La méthode d'Euler progressive devient instable à partir d'un pas "
         "qui vaut %f\n",
         2.0 / pvals[0]);
  double h = (2.0 / pvals[0]) - (2 / (pvals[0] * 100));
  printf("Configurons un pas de %f\n\n", h);

  for (int i = 0; i < n; i++) {
    U[i] = T0;
  }

  FILE *fp = NULL; // Ouvrir le fichier pour écrire les données
  fp = popen("gnuplot -persist", "w");
  if (fp) {
    fprintf(fp, "unset key\n");
    fprintf(fp, "set view map\n");
    fprintf(fp, "set cbrange [0:20]\n");

    for (int k = 0; k < (2000.0 / h); k++) {
      if (k % m == 0) {
        fprintf(fp, "splot '-' with pm3d\n");
        for (int i = 0; i < ne; i++) {
          if ((3.0 <= datax[i]) && (datax[i] <= 6.0) && (3.0 <= datay[i]) &&
              (datay[i] <= 6.0)) { /*On vérifie si on se trouve dans le trou ou
                                      sur le bord*/
            error++;
            fprintf(fp, "%f %f %f\n", datay[i], datax[i], 0.0);
            continue;
          } else if ((datay[i] == 0.0) || (datay[i] == 10.0) ||
                     (datax[i] == 0.0) ||
                     (datax[i] == 8.0)) { /*On vérifie si on se trouve sur le
                                             bord de la grille*/
            fprintf(fp, "%f %f %f\n", datay[i], datax[i], 0.0);
            error++;
          } else {
            fprintf(fp, "%f %f %f\n", datay[i], datax[i], U[i - error]);
          }

          if ((i + 1) % (nx + 2) == 0) {
            fprintf(fp, "\n");
          }
        }
        fprintf(fp, "e\n");

        fflush(fp); // Forcer l'actualisation du graphique
        if (m <= 2)
          usleep(500000);
        else
          usleep(100000);
      }
      matvec(U, result, &nev, n, ia, ja,
             a); /*On calcul le résultat du prorduit matrice vecteur entre A et
                    U(t)*/

      for (int j = 0; j < n; j++) {
        result[j] = -D * result[j];
        U[j] = U[j] + (h * result[j]); /*On recalcul U pour pouvoir l'afficher à
                                          l'itération suivante*/
      }
      error = 0.0;
    }

    pclose(fp);
  }
  /*Calcul de la valeur propre et du vecteur propre par EVSL*/

  double *evsl_eval = malloc(nev * sizeof(double));
  double *evsl_evecs = malloc(nev * n * sizeof(double));

  int ierr = evsl(m, n, ia, ja, a, nev, evsl_eval, evsl_evecs);
  if (ierr != 0) {
    printf("Erreur lors du calcul des valeurs propres et vecteurs propres par "
           "EVSL\n");
  } else {
    printf("Valeur propre minimale calculée par EVSL: %f\n", evsl_eval[0]);
  }

  /*libérer la mémoire */
  free(ia);
  free(ja);
  free(a);
  free(evals);
  free(evecs);
  free(datax);
  free(datay);
  free(U);
  free(result);
  free(pvals);
  free(pvecs);
  free(Da);
  free(evsl_eval);
  free(evsl_evecs);
  return 0;
}
