#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define BLEN 10
#define DATA_TAG 99

int main (int argc, char **argv)
{

    int i, j;
    int *data, *buffer;
    int my_rank, num_procs;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    /* MASTER PART */
    if(my_rank == 0){

      /* allocate memory */
      data = (int*) malloc (BLEN * num_procs * sizeof(int));
      buffer = (int*) malloc (BLEN * num_procs * sizeof(int));

      /* initialize array */
      for(i=0; i<BLEN*num_procs; ++i){
	data[i] = i;
	/*printf(" %d ",data[i]);*/
      }

      /* prepare buffer for each processor and send */
      for(i=0; i<num_procs; ++i){
	for(j=0; j<BLEN; ++j){
	  buffer[i*BLEN+j] = data[j*num_procs+i]; 
	}
	if(i == 0){
	  printf("Process %d recieved: ",my_rank);
	  for(j=0; j<BLEN; j++)
	    printf(" %d ",buffer[j]);
	  printf("\n");
	}
	else
          MPI_Send(&buffer[i*BLEN], BLEN, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD);
      }
    }

    /* WORKER PART */
    else{
      buffer = (int*) malloc (BLEN * sizeof(int));
      MPI_Recv(buffer, BLEN, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, &status);
      printf("Process %d recieved: ",my_rank);
      for(j=0; j<BLEN; j++)
	printf(" %d ",buffer[j]);
      printf("\n");
    }

    /* CLEAN UP */
    if (my_rank == 0)
      free(data);
    free(buffer);
    MPI_Finalize();

    return EXIT_SUCCESS;
}
