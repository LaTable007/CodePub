#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define P_SEP printf("====================================\n")
int prob(int m, int *n, int **ia, int **ja, double **a)
/*
   But
   ===
   Génère la matrice n x n qui correspond à la discrétisation sur une grille 
   cartésienne régulière m x m de l'opérateur de Laplace à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )        sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec la fonction u qui satisfait les conditions aux limites de Dirichlet
         
         u = 0  sur (0,y), (1,y), (x,0) et (x,1), avec 0 <= x,y <= 1 .
  
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
{   
  int  nnz, ix, iy, nx, ind = 0, ny;
  double invh2, lx = 2.0, ly = 5.0/2, d, h; /* longueur du côté le plus petit du rectangle (en m) */

  if (lx <= ly){ // si le rectangle est plus long en y on fait en sorte que le pas de discrétisation est définie sur x
    nx = m - 2; //nombre d'élément le long de x
    d = ly - lx; // la différence entre la longueur et largeur
    h = lx / (m-1);// le pas de discrétisation
    ny = (d / h) + nx;//nombre d'élément selon y
    invh2 = (m-1) * (m-1) / (lx * lx);//inverse au carré du pas de dicrétisation
  }

  else{//même résonnement qu'avant mais pour un rectangle dont le plus grand côté est en x
    ny = m - 2;
    d = lx - ly;
    h = ly / (m-1);
    nx = (d / h) + ny;
    invh2 = (m-1) * (m-1) / (ly * ly);
  }

  printf("le pas = %f\n", h);
  printf("élément x = %d\n", nx);
  printf("élément y = %d\n", ny);
  printf("invh2 = %f\n", invh2);

  double tx1 = 3.0/4, tx2 = 6.0/4, ty1 = 3.0/4, ty2 = 6.0/4;
  int tind, dx, find, dy;

  if ((tx1 != 0.0) && (ty1 != 0.0) && (tx2 != lx)){
    printf("cas 1\n");
    tind = (tx1 / h) + (nx * ((ty1 / h) - 1)) - 1;
    dx = (tx2 / h) + (nx * ((ty1 / h) - 1)) - tind;
    dy = (((ty2 - ty1) / h) * nx);
    find = (tx2 / h) + (nx * ((ty2 / h) -1)) - 1;  
  }
/*
  else if ((tx1 != 0.0) && (ty1 == 0.0)) {
    printf("cas 2\n");
    tind = tx1 / h;
    dx = tx2 / h - tind;
  }

  else if ((tx1 == 0.0) && (ty1 != 0.0)) {
    printf("cas 3\n");
    tind = 1 + (nx * ((ty1 / h) - 1));
    dx = (tx2 / h) + (nx * ((ty1 / h) - 1)) - tind;
  }

  else if ((tx1 == 0.0) && (ty1 == 0.0)) {
    printf("cas 4\n");
    tind = 0.0;
    dx = (tx2 / h) + (nx * (ty1 / h)) - 1.0;
  }

  else if (tx2 == lx){
    printf("cas 5\n");
    tind = (tx1 / h) + (nx * ((ty1 / h) - 1));
    dx = (tx2 / h) + (nx * ((ty1 / h) - 1)) - tind - 1;
  }
*/
  printf("indice trou = %d\n", tind);

  printf("indice fin = %d\n", find);

  printf("largeur trou = %d\n", dx);

  printf("longueur trou = %d\n", dy);

  *n  = nx * ny - (((tx2 - tx1 + h) / h) * (ty2 - ty1 + h) / h); /* nombre d'inconnues */
  nnz = ((5 * nx * ny) - (4.0 * (ny + nx) / 2)) - (6.0 * (tx2 - tx1 + h) * (ty2 - ty1 + h) / (h * h));
  /* nombre d'éléments non nuls a retiré calculé 64 normalement mais on retire 96 pour arriver à 187  */
  
  printf("nnz = %d\n", nnz);
  printf("le nombre d'élément = %d\n", *n);
  
  /* allocation des tableaux */

  *ia  = malloc((*n + 1) * sizeof(int));
  *ja  = malloc(nnz * sizeof(int));
  *a   = malloc(nnz * sizeof(double));
  
  /* allocation réussite? */

  if (*ia == NULL || *ja == NULL || *a == NULL ) {
    printf("\n ERREUR : pas assez de mémoire pour générer la matrice\n\n");
    return 1;
  }

  /* partie principale : remplissage de la matrice */

  nnz = 0;
  int error = 0;
  int vind = 0;

  for (iy = 0; iy < ny; iy++) {
    for (ix = 0; ix < nx; ix++) { 
      /* numéro de l'équation */
      ind = ix  + nx * iy;

      //printf("%d ", ind); 
      
      if ((tind <= ind) && ((ind - tind) % nx <= dx - 1) && (ind <= find)) {// L'indice ind est dans le trou carré
        //printf("L'indice %d se trouve dans le trou carré.\n", ind);
        error++;
        continue;
      }
      /*
      else {// L'indice ind n'est pas dans le trou carré
        printf("L'indice %d n'est pas dans le trou carré.\n", ind);
      }
      */

      //printf("%d ", error);

      vind = ind - error;

      //printf("%d ", vind);

      //printf("%d", vind);

      /* marquer le début de la ligne suivante dans le tableau 'ia' */
      (*ia)[vind] = nnz;       
      //printf("%d\n", (*ia)[vind]);
        
      /* remplissage de la ligne : voisin sud */
      if (iy > 0){
        if ((ind < find - dx + nx + 1) || (find + nx < ind))  {
          (*a)[nnz] = -invh2; /* pour D=1 */
          if((ind < tind) || (find +nx - dx < ind )){
            (*ja)[nnz] = vind - nx;
          }
          else{
            (*ja)[nnz] = vind - nx + dx;
          }
          //printf("Sud %d / %f / %d\n", vind, (*a)[nnz], (*ja)[nnz]);
          nnz++;
        }
      }

      /* remplissage de la ligne : voisin ouest */
      if (ix > 0){
        if((tind + dx <= ind) && ((ind - tind - dx) % nx == 0) && (ind <= find + 1)){
        }
        else {
          (*a)[nnz] = -invh2; /* pour D=1 */
          (*ja)[nnz] = vind - 1;
          //printf("Ouest %d / %f / %d\n", vind, (*a)[nnz], (*ja)[nnz]);
          nnz++;
        }
      }
      
      /* remplissage de la ligne : élément diagonal */
      (*a)[nnz] = 4.0*invh2; /* pour D=1 */
      (*ja)[nnz] = vind;
      //printf("Diagonal %d / %f / %d\n", vind, (*a)[nnz], (*ja)[nnz]);
      nnz++;

      /* remplissage de la ligne : voisin est */
      if (ix < nx - 1 ){
        if((tind - 1<= ind) && ((ind - tind + 1) % nx == 0) && (ind <= find - dx)){
        }
        else{
          (*a)[nnz] = -invh2; /* pour D=1 */
          (*ja)[nnz] = vind + 1;
          //printf("Est %d / %f / %d\n", vind, (*a)[nnz], (*ja)[nnz]);
          nnz++;
        }
      }

      /* remplissage de la ligne : voisin nord */
      if (iy < ny - 1){
        if((ind < tind - nx) || (ind > tind + dx - 1 - nx)) {
          (*a)[nnz] = -invh2; /* pour D=1 */
          if((ind < tind + dx - nx) || (find - dx < ind)){
            (*ja)[nnz] = vind + nx;
          }
          else{
            (*ja)[nnz] = vind + nx - dx;
          }
          //printf("Nord %d / %f / %d\n", vind, (*a)[nnz], (*ja)[nnz]);
          nnz++; 
        } 
      }
      //printf("\n");
    }
  } 

  /* dernier élément du tableau 'ia' */
  (*ia)[vind + 1] = nnz;

  P_SEP;
  int f = 0;
  for(;f < nnz;f++){
   printf("%f\n", (*a)[f]);
  }
  P_SEP;
  int i = 0;
  for(;i < *n + 1; i++) {
    printf("%i\n", (*ia)[i]);
  }
  P_SEP;
  int j = 0;
  for(;j < nnz; j++) {
    printf("%d\n", (*ja)[j]);
  }

  /* retour habituel de fonction */
  return 0;
}
