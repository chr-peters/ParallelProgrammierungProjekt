/*  Serial Program */
/*  Approximate pi with the n-point rectangle quadrature rule */
/*  applied to the integral from 0 to 1 of 4 / (1+x**2).c */

#include <stdio.h>

int main(int argc, char *argv[])
{
  FILE *in;
  
  double pi  = 0.0; /* The calculated result */
  int    n   = 0;   /* Number of points of integration */
  double h;         /* Reciprocal of n, interval size */
  double x;         /* Midpoint of each rectangle's interval */
  double f;         /* Function to integrate */
  double sum;       /* Area of rectangles */

  int    i; /* do loop index */
  
  in = fopen("input_pi.dat", "r");
  
  if (in == NULL)
    printf("ERROR: Cannot open file input_pi.dat!\n");
  else{
    fscanf(in,"%d",&n);
    fclose(in);
  }
  
  if(n > 0) {
    
    h   = 1.0 / n; /* Calculate interval size */
    sum = 0.0;     /* Initialize sum */
    
    /* Calculate partial sums */
    for(i = 1; i <= n; i++) {
      x   = (i - 0.5) * h;
      f   = 4.0 / (1.0 + x * x);
      sum = sum + f;
    }
    
    pi = h * sum;
    
    printf("Value of pi is: %20.16lf\n", pi);
  }
  else {
    printf("ERROR: %d is not a valid value for n. n must be > 0\n",n);
  }
  return 0;
}
