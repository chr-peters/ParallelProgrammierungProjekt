#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define BLOCKSIZE 10000000
#define MEASUREMENTS 20

double esecond(void) {

  struct timeval tp;
  struct timezone tzp;
  
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + (tp.tv_usec * 1e-6);
}

int main() {
  char * sched_types[]={"serial", "static", "dynamic", "guided", "auto"};
  static float a[BLOCKSIZE];
  double s, t;
  double start, end;
  int nthreads, i;
  int sched_type;
  int chunksize;
  int sweep;

  start = esecond();

  for (sweep=0; sweep < MEASUREMENTS; sweep++)
  {
    /* TODO: set variables */
    nthreads = ...
    sched_type = ...
    chunksize = ...

    /* TODO: parallelize loop with OpenMP */  
    for(i=0; i<BLOCKSIZE; i++) {
      a[i] = i*17;
    }

    s=0.0;

    /* TODO: parallelize loop with OpenMP */  
    for(i=0; i<BLOCKSIZE; i++) {
      t = a[i];
      s = s + sqrt(t);
    }
  } /* end for */

  end = esecond();

  printf("... %10.3f seconds on %d pe(s) (type: %s, chunksize: %d)\n",
         (end - start)/MEASUREMENTS, nthreads, sched_types[sched_type], chunksize);

  return(0);
}

