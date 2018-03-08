/* ---------------------------------------------------------------- */
/* Serielles Program zur Bestimmung des Vorkommens einer Zahl in    */
/* einem Integer-Vektor                                             */
/* ---------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>

int main (int argc, char **argv){

    FILE *in;
  
    int i, k;
    int count = 0;
    int number;

    int *iVector;

    /* initialize */
    in = fopen("input_numarray.dat", "r");
    
    if (in == NULL)
      printf("ERROR: Cannot open file input_numarray.dat!\n");
    else{
      fscanf(in,"%d",&k);
      fscanf(in,"%d",&number);
      fclose(in);
    }

    
    iVector = (int*) malloc (k * sizeof(int));

    /* fill array with random numbers between 0 and 9 */
    printf("\nk=%d  -->  search number = %d\n",k,number);

    printf("\nArray:\n");
    for(i=0; i<k; i++){
        iVector[i] = rand()%10;
        printf(" %i", iVector[i]);
        if((i+1)%20 == 0) printf("\n");
    }

    printf("\n");

    /* count appearance of search number */
    for(i=0; i<k; i++) {
        if(iVector[i] == number) count ++;
    }

    printf("\nNumber %d was found %d times in the array\n\n", number, count);

    free(iVector);

    return 0;
}
