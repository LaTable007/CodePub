#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define P_SEP printf("====================================\n")
int prob(int m, int *n, int **ia, int **ja, double **a, double **datax,
         double **datay, int *ne, int *nx)
/*
   But
   ===
   Génère la matrice n x n qui correspond à la discrétisation sur une grille
   cartésienne régulière nx x ny de l'opérateur de Laplace à deux dimension.

   Dans notre cas le m sera définie à partir d'un carré de longueur unité
   donc les longueur fournie dans l'énoncé du projet ont été multiplié par 4
  afin de satisfaire cette condition. Il faura faire de même si le lecteur
  souhaite changer les dimensions de la grille.

  La numérotation des inconnues est lexicographique, la direction x étant
  parcourue avant celle de y. La matrice est retournée dans le format CRS
  qui est défini par le scalaire 'n' et les trois tableaux 'ia, 'ja' et 'a'.

  Arguments
  =========
  m (input)   - nombre de points par direction dans la grille
  n  (output) - pointeur vers le nombre d'inconnues dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A

  Sortie
  ======
  0 - exécution avec succès
  1 - erreurs
*/
{ /*calcul des valeurs du problème*/
  int nnz, ix, iy, ind = 0, ny;
  double invh2, lx = 8.0, ly = 10.0, d, h;

  *nx = (lx * m) - 1; /*nombre d'élément dans la direction x*/
  d = ly - lx;
  h = 1.0 / (m);
  ny = (d / h) + *nx; /*nombre d'élément dans la direction y*/
  invh2 = 16.0 / (h * h);

  double tx1 = 3.0, tx2 = 6.0, ty1 = 3.0, ty2 = 6.0;
  int tind, dx, find, dy;

  tind = (tx1 / h) + (*nx * ((ty1 / h) - 1)) - 1; /*indice du début du trou*/
  dx = (tx2 / h) + (*nx * ((ty1 / h) - 1)) -
       tind;                  /*nombre d'élément dans la direction x du trou*/
  dy = ((ty2 - ty1) / h) + 1; /*nombre d'élément dans la direction y du trou*/
  find = (tx2 / h) + (*nx * ((ty2 / h) - 1)) - 1; /*indice de fin du trou*/

  *ne = (*nx + 2) * (ny + 2); /*nombre d'élément dans la grille en contant le
                                 trou (sera utilisé pour la question 3 et 4)*/
  *n = *nx * ny - (((tx2 - tx1 + h) / h) * (ty2 - ty1 + h) / h);
  nnz = (5 * *nx * ny) - (2 * (*nx + ny)) - (5 * dx * dy) - (2 * (dy + dx));

  *ia = malloc((*n + 1) * sizeof(int));
  *ja = malloc(nnz * sizeof(int));
  *a = malloc(nnz * sizeof(double));

  if (*ia == NULL || *ja == NULL || *a == NULL) {
    printf("\n ERREUR : pas assez de mémoire pour générer la matrice\n\n");
    return 1;
  }

  nnz = 0;
  int error = 0;
  int vind = 0;

  /*calcul des coordonnées x et y pour la question 3 et 4*/
  *datax = malloc((*ne) * sizeof(double));
  *datay = malloc((*ne) * sizeof(double));

  if (*datax == NULL || *datay == NULL) {
    printf("\n ERREUR : pas assez de mémoire pour générer les coordonnées\n\n");
    return 1;
  }

  for (iy = 0; iy < (ny + 2); iy++) {
    for (ix = 0; ix < (*nx + 2); ix++) {
      int dataind = ix + ((*nx + 2) * iy);
      (*datax)[dataind] = h * (ix);
      (*datay)[dataind] = h * (iy);
    }
  }

  /*remplissage de la matrice dans le format CSR*/

  for (iy = 0; iy < ny; iy++) {
    for (ix = 0; ix < *nx; ix++) {
      ind = ix + *nx * iy;

      if ((tind <= ind) && ((ind - tind) % *nx <= dx - 1) && (ind <= find)) {
        error++; /*On utlise error afin de ne pas prendre en compte les éléments
                    qui se trouvent dans le trou dans la matrice A*/
        continue;
      }

      vind = ind - error; /*Représente l'indice de la matrice A*/

      (*ia)[vind] = nnz;
      /*Remplissage de la case voisin Sud*/
      if (iy > 0) {
        if ((ind < find - dx + *nx + 1) ||
            (find + *nx <
             ind)) { /*On vérifie qu'on ne se trouve pas au dessus du trou*/
          (*a)[nnz] = -invh2;
          if ((ind < tind) ||
              (find + *nx - dx <
               ind)) { /*On vérifie qu'on est pas encore passé dans le trou*/
            (*ja)[nnz] = vind - *nx;
          } else {
            (*ja)[nnz] = vind - *nx + dx;
          }
          nnz++;
        }
      }
      /*Remplissage de la ligne voisin Ouest*/
      if (ix > 0) {
        if ((tind + dx <= ind) && ((ind - tind - dx) % *nx == 0) &&
            (ind <=
             find + 1)) { /*On vérifie qu'on ne se trouve pas à droite du trou*/
        } else {
          (*a)[nnz] = -invh2;
          (*ja)[nnz] = vind - 1;
          nnz++;
        }
      }
      /*Remplissage de la ligne voisin Nord*/
      (*a)[nnz] = 4.0 * invh2;
      (*ja)[nnz] = vind;
      nnz++;
      /*Remplissage de la ligne voisin Est*/
      if (ix < *nx - 1) {
        if ((tind - 1 <= ind) && ((ind - tind + 1) % *nx == 0) &&
            (ind <=
             find -
                 dx)) { /*On vérifie qu'on ne se trouve pas à gauche du trou*/
        } else {
          (*a)[nnz] = -invh2;
          (*ja)[nnz] = vind + 1;
          nnz++;
        }
      }
      /*Remplissage de la ligne voisin Sud*/
      if (iy < ny - 1) {
        if ((ind < tind - *nx) ||
            (ind > tind + dx - 1 - *nx)) { /*On vérifie qu'on ne se trouve pas
                                              en dessous du trou*/
          (*a)[nnz] = -invh2;
          if ((ind < tind + dx - *nx) ||
              (find - dx <
               ind)) { /*On vérifie qu'on n'est pas encore passé dans le trou*/
            (*ja)[nnz] = vind + *nx;
          } else {
            (*ja)[nnz] = vind + *nx - dx;
          }
          nnz++;
        }
      }
    }
  }

  (*ia)[vind + 1] = nnz;
  return 0;
}
