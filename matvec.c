#include <stdlib.h>
#include <stdio.h>

/* variables statiques -- accessibles  */

void matvec(void *vx, void *vy, int *blockSize, int n, int *ia, int *ja, double *a)
/*
   But
   ===
   Calcule le produit matrice-vecteur
                              vy = A*vx 
   pour le solveur aux valeurs propres PRIMME. La matrice A doit être
   "stoquée" au préalable dans les variables statiques 'n', 'ia', 'ja' et 'a'
   en utilisant le format CSR (Compressed Sparse Rows). Par "stoquer" 
   on veut dire ici stoquer la valeur de 'n' et les pointeurs vers les 
   tableaux 'ia', 'ja' et 'a'.

   Arguments
   =========
   vx        (input) - vecteur(s) d'entrée
   vy       (output) - vecteur(s) de produit A*vx
   blockSize (input) - nombre de vecteurs d'entrée
   primme    (input) - paramètres fournis par primme pour optimiser le calcul
                       (pas utilisé) 
*/
{
    int i, j, b;
    double *x = vx, *y = vy;

    for(b = 0; b < (*blockSize)*n; b+=n)
        for(i = 0; i < n; i++){
            y[b+i] = 0;
            for (j = ia[i]; j < ia[i + 1]; j++)
                y[b+i] += a[j] * x[b+ja[j]];
        }
} 
