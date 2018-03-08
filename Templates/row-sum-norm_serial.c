/******************************************************************************/
/*  Serial programm for calculating the row-sum norm of a (n x n) matrix      */
/******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAG 101

#define N 16  /* matrix dimesion */

void init_matrix(int *matA);

void calc_row_sums(int *matA, int dim, int *sums);

void check_result(int *matA);

/******************************************************************************/
/*                              M A I N                                       */
/******************************************************************************/

int main(int argc, char *argv[]) 
{

  int *matA;                 /* global matrix */
  int *sums;                  /* global row sums */
  int norm;
  
  int i, j;

  /* initialize matrix */
  matA = (int *) malloc(N * N * sizeof(int));
  init_matrix(matA);

  /* calculate row sums */
  sums = (int *) malloc(N * sizeof(int));
  calc_row_sums(matA, N, sums);

  /* determine maximum of row sums */
  norm = 0;
  for (i=0; i<N; ++i){
    if (sums[i] > norm)
      norm = sums[i];
  }
  printf("row-sum norm: %d\n", norm);

  /*  check correctness */
  check_result(matA);

  return 0;
}


/******************************************************************************/
/*  init matrices                                                             */
/******************************************************************************/
void init_matrix(int *matA){

  int i, j, k;
  int num_test_prints;

  for (i=0; i<N; i++) {
    for (j=0;j<N; j++) { 
      matA[i*N+j]= -10 + rand()%20; 
    }
  }
  /* test prints */
  for (i=0;i<N; i++) { 
    for (j=0;j<N; j++) { 
      printf("%3d",matA[i*N+j]);
    }
    printf("\n");
  }
}


/******************************************************************************/
/*  calc_row_sums                                                             */
/******************************************************************************/
void calc_row_sums(int *matA, int dim, int *sums){

  int i, j;

  for (i=0; i<dim; i++) {
    sums[i] = 0;
    for (j=0;j<dim; j++) {
      sums[i] += abs(matA[i*dim+j]);
    }
  }
}

/******************************************************************************/
/*  init matrices                                                             */
/******************************************************************************/
void check_result(int *matA){

  int i, j;
  int row_sum;
  int norm=0;

  for (i=0; i<N; i++) {
    row_sum = 0;
    for (j=0;j<N; j++) {
      row_sum += abs(matA[i*N+j]);
    }
    if (row_sum > norm)
      norm = row_sum;
  }

  printf("serial row-sum norm: %d\n",norm);
}
