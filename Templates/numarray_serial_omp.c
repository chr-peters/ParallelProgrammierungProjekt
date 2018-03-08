/* ------------------------------------------------------ */
/* Serielles Program welches das Vorkommen der 0 in einem */       
/* Integer-Vektor zählt                                   */
/* ------------------------------------------------------ */

#include <stdlib.h>
#include <stdio.h>

int main (int argc, char **argv){

    int i, k;
    int count = 0;

    int *iVector;

    printf("Enter array dimension: ");
    scanf("%d",&k);

    iVector = (int*) malloc (k * sizeof(int));

    /* fill array with random numbers */
    printf("\nk=%d\n",k);

    printf("\nArray:\n");
    for(i=0; i<k; i++){
        iVector[i] = rand()%10;
        printf(" %i", iVector[i]);
        if((i+1)%20 == 0) printf("\n");
    }

    printf("\n");

    /* count zeros */
    for(i=0; i<k; i++) {
      if(iVector[i] == 0) 
	count ++;
    }

    printf("\nNumber of zeros: %d \n\n", count);

    free(iVector);

    return 0;
}
