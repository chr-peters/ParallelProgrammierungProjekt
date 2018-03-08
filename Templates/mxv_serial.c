/**************************************/
/*   Matrix-Vektor-Multiplikation     */          
/**************************************/

#include <stdio.h>

#define SIZE 10


int main (void){

  float A[SIZE][SIZE], b[SIZE], c[SIZE], total;
  int i, j;
  
  /* Initializations */
  total = 0.0;
  for (i=0; i < SIZE; i++)
    {
      for (j=0; j < SIZE; j++)
	A[i][j] = (j+1) * 1.0;
      b[i] = 1.0 * (i+1);
      c[i] = 0.0;
    }
  printf("\nStarting values of matrix A and vector b:\n");
  for (i=0; i < SIZE; i++)
    {
      printf("  A[%d]= ",i);
      for (j=0; j < SIZE; j++)
	printf("%.1f ",A[i][j]);
      printf("  b[%d]= %.1f\n",i,b[i]);
    }
  printf("\nIntermediate results:\n");
  
  for (i=0; i < SIZE; i++)
    {
      for (j=0; j < SIZE; j++)
	c[i] += (A[i][j] * b[j]);
      
      total = total + c[i];
/*       printf("Row %d\t c[%d]=%.2f\t",i,i,c[i]); */
/*       printf("Running total= %.2f\n",total); */
    }   /* end of parallel i loop */
  
  printf("\nMatrix-vector total - sum of all c[] = %.2f\n\n",total);
  return 0;
}
