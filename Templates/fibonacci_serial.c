/**********************************************************/
/* Serielles Programm zur Berechnung einer Fibonacci-Zahl */
/**********************************************************/

#include <stdio.h>

int fib(int n){

  int i,j;
  if (n<2)
    return n;
  else{
    i=fib(n-1);
    j=fib(n-2);
    return i+j;
  }
}

int main(void){

  int num_fib;

  printf("Which fibonacci number do you want to compute?\n");
  scanf("%d",&num_fib);

  printf("\nThe %d'th Fibonacci number is: %d\n", num_fib, fib(num_fib));

  return 0;
}
