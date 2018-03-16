/*
 * mandelseq.c
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <ppmwrite.h>
#include <mpi.h>
#include <math.h>

#define WORK_TAG 69
#define DIE_TAG 88

double esecond(void) {

  struct timeval tp;
  struct timezone tzp;

  gettimeofday(&tp, &tzp);
  return tp.tv_sec + (tp.tv_usec * 1e-6);
}

void usage(char *name) {
  fprintf(stderr, "Usage: %s options\n\nwith the following optional options (default values in parathesis):\n\n",name);

  fprintf(stderr, "  [-x <x0> <x1> <y0> <y1>]  coordinates of initial area (-1.5 0.5 -1.0 1.0)\n");
  fprintf(stderr, "  [-w <width>]              image width in pixels (256)\n");
  fprintf(stderr, "  [-h <height>]             image height in pixels (256)\n");
  fprintf(stderr, "  [-i <maxiter>]            max. number of iterations per pixel (256)\n");
  fprintf(stderr, "  [-t <type>]               0=stride, 1=stripe, 2=blockmaster\n");
  fprintf(stderr, "  [-v]                      verbose (off)\n\n");
  exit(1);
}

void calc(int *iterations, int width, int height, int xstart, int ystart, double xmin, double xmax, double ymin, double ymax, int maxiter);

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  /* values for MPI */
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* default values for command line parameters */
  double xmin = -1.5;  /* coordinates of rectangle */
  double xmax =  0.5;
  double ymin = -1.0;
  double ymax =  1.0;
  int width   = 256;   /* size of rectangle in pixels */
  int height  = 256;
  int maxiter = 256;   /* max. number of iterations */
  int verbose = 0;     /* per default only print error messages */
  int type = 0;        /* per default only print error messages */
  int *iterations,*recvbuffer;
  int ix,iy;
  int i;

  double calctime=0.0, runtime=0.0, mpitime=0.0, waittime=0.0, iotime=0.0;
  char   filename[1024];

  // measure runtime
  MPI_Barrier(MPI_COMM_WORLD);
  double runtime_start = esecond();

  ppminitsmooth(1);

  /* parse command line */
  i=1;
  while( i < argc ) {
    if( argv[i][0] == '-' ) {
      switch( argv[i][1] ) {
      case 'x':
        xmin = atof(argv[++i]);
        xmax = atof(argv[++i]);
        ymin = atof(argv[++i]);
        ymax = atof(argv[++i]);
        break;
      case 'i':
        maxiter = atoi(argv[++i]);
        break;
      case 'w':
        width = atoi(argv[++i]);
        break;
      case 'h':
        height = atoi(argv[++i]);
        break;
      case 't':
        type = atoi(argv[++i]);
        break;
      case 'v':
        verbose++;
        break;
      default:
        usage(argv[0]);
      }
    } else {
      usage(argv[0]);
    }
    i++;
  }

  /* initialize arrays */
  if (rank == 0) {
    recvbuffer = malloc(width*height*sizeof(int));
    for (ix=0; ix<width; ++ix) {
      for (iy=0; iy<height; ++iy) {
	       recvbuffer[ix*height+iy] = 0;
      }
    }
  }

  if (type == 0) {
    /**
     * Initialize memory.
     */
    iterations = malloc(width*height*sizeof(int));

    for (ix=0; ix<width; ++ix) {
      for (iy=0; iy<height; ++iy) {
	       iterations[ix*height+iy] = 0;
      }
    }

    // starting measurement of calctime
    double calcstart = esecond();

    /**
     * Perform calculations.
     */
    double rowHeight = (ymax-ymin)/height;
    int curRow;
    for (curRow=rank; curRow<height; curRow+=size){
      calc(iterations, width, 1, 0, curRow, xmin, xmax, ymin + curRow*rowHeight, ymin + curRow*rowHeight+rowHeight, maxiter);
    }

    // measure calctime
    calctime = esecond() - calcstart;

    // measure wait time
    double waitstart = esecond();
    MPI_Barrier(MPI_COMM_WORLD);
    waittime = esecond() - waitstart;

    /**
     * Get the results.
     */
    // starting measurement of mpitime
    double mpistart = esecond();
    MPI_Reduce(iterations, recvbuffer, width*height, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    mpitime = esecond() - mpistart;
  }

  else if (type == 1) {
    // calculate number of rows per process
    int rowCount = ceil((float)height / size);
    rowCount = fmin(rowCount, height - rank * rowCount);

    // initialize the iterations field
    iterations = malloc(width*rowCount*sizeof(int));

    // calculate the y range of the block
    double yRange = ((double)rowCount/height)*(ymax-ymin);

    int defaultRowCount = ceil((float)height / size);
    double defaultYRange = ((double)defaultRowCount/height)*(ymax-ymin);

    // starting measurement of calctime
    double calcstart = esecond();

    // call the calc method
    calc(iterations, width, rowCount, 0, 0, xmin, xmax, ymin + rank*defaultYRange, ymin + rank*defaultYRange + yRange, maxiter);

    // measure calctime
    calctime = esecond() - calcstart;

    // measure wait time
    double waitstart = esecond();
    MPI_Barrier(MPI_COMM_WORLD);
    waittime = esecond() - waitstart;

    /**
     * Collect the results
     */
    int *recvcounts;
    int *displs;
    if (rank == 0) {
      // only initialize the arrays on the root process
      recvcounts = malloc(size * sizeof(int));
      displs = malloc(size * sizeof(int));

      // fill the arrays
      //int i;
      // the default row count of each process
      int defaultRowCount = ceil((float)height / size);
      for (i=0; i<size; i++) {
      	if (i==size-1) {
      	  int tmpRowCount = fmin(defaultRowCount, height - i * defaultRowCount);
      	  recvcounts[i] = tmpRowCount * width;
      	  displs[i] = defaultRowCount * width * i;
      	} else {
      	  recvcounts[i] = defaultRowCount * width;
      	  displs[i] = defaultRowCount * width * i;
      	}
      }
    }

    // measure mpi time
    double mpistart = esecond();
    MPI_Gatherv(iterations, rowCount*width, MPI_INT, recvbuffer, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    mpitime = esecond() - mpistart;
  }

  else if ( type == 2 ) {
      int blocksize;
      /*
       *  Find optimal dimension for the blocks, starting at 16x16
       */
      for(i = 0; i < 15; ++i){
        if((height % (16-i) == 0) && (width % (16-i) == 0)){
          blocksize = 16-i;
          break;
        }
      }
      //if(verbose && rank == 0) printf("Blocksize: %d\n", blocksize);
      /*
       *  Define vector for optimal blocksize
       */
      MPI_Datatype vector;
      MPI_Type_vector(blocksize, blocksize, width, MPI_INT, &vector);
      MPI_Type_commit(&vector);
      MPI_Status status;

      /*
       * Master
       */
      if (rank == 0) {
          int send_block_number = 0;
          int receive_block_number;
          int num_tasks = width * height / (blocksize * blocksize);
          // Send initial tasks for slaves
          for (i = 1; i < size; ++i) {
              if (send_block_number < num_tasks) {
                  MPI_Send(&send_block_number, 1, MPI_INT, i, WORK_TAG, MPI_COMM_WORLD);
                  send_block_number++;
              }
          }

          // Receive data from slaves and send the remaining tasks individually
          int offset_x;
          int offset_y;
          while (send_block_number < num_tasks) {
              MPI_Recv(&receive_block_number, 1, MPI_INT, MPI_ANY_SOURCE, WORK_TAG, MPI_COMM_WORLD, &status);
              offset_x = (receive_block_number * blocksize) % width;
              offset_y = (receive_block_number * blocksize) / width * width * blocksize;
              MPI_Recv(recvbuffer + offset_x + offset_y, 1, vector, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD, &status);
              MPI_Send(&send_block_number, 1, MPI_INT, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD);
              send_block_number++;
          }

          // Receive remaining data after all tasks are sent
          for (i = 1; i < size; ++i) {
              MPI_Recv(&receive_block_number, 1, MPI_INT, MPI_ANY_SOURCE, WORK_TAG, MPI_COMM_WORLD, &status);
              offset_x = (receive_block_number * blocksize) % width;
              offset_y = (receive_block_number * blocksize) / width * width * blocksize;
              MPI_Recv(recvbuffer + offset_x + offset_y, 1, vector, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD, &status);
          }

          // Kill the slaves
          int dummy = 0;
          for (i = 1; i < size; ++i) {
              MPI_Send(&dummy, 1, MPI_INT, i, DIE_TAG, MPI_COMM_WORLD);
          }

      /*
       * Slave
       */
      } else {
          int block_number;
          iterations = malloc(sizeof(int) * blocksize * blocksize);

          while (1) {
              MPI_Recv(&block_number, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              if (status.MPI_TAG == DIE_TAG){
                break;
              }
              //Calculate values for calc method
              double xrange = xmax - xmin;
              double xlower = xmin + ((block_number * blocksize) % width) * xrange/width;
              double xupper = xlower + (double) blocksize/width * xrange;
              double yrange = ymax - ymin;
              //printf("ymax: %lf ymin: %lf\n", ymax, ymin);
              double ylower = ymin + (double)((block_number * blocksize)/width * blocksize)/height * yrange;
              double yupper = ylower + (double) blocksize/height * yrange;
              /*printf("Blocksize: %d Height: %d yrange: %lf\n", blocksize, height, yrange);
              printf("xrange: %lf yrange: %lf\n", xrange, yrange);*/
              //printf("x -> [%lf %lf] y -> [%lf %lf]\n", xlower, xupper, ylower, yupper);
              calc(iterations, blocksize, blocksize, 0, 0, xlower, xupper, ylower, yupper, maxiter);
              MPI_Send(&block_number, 1, MPI_INT, 0, WORK_TAG, MPI_COMM_WORLD);
              MPI_Send(iterations, blocksize*blocksize, MPI_INT, 0, WORK_TAG, MPI_COMM_WORLD);
          }
      }
      MPI_Type_free(&vector);
  }

  // measure time of whole program
  MPI_Barrier(MPI_COMM_WORLD);
  runtime = esecond() - runtime_start;

  // output the picture
  if(rank == 0){//master
    // measure iotime
    double iostart = esecond();
    ppmwrite(recvbuffer,width,height,0,maxiter,"mandelcol.ppm");
    iotime = esecond() - iostart;
  }

  // print the time data in a csv line
  if (verbose) {
    printf("%d,%d,%d,%f,%f,%f,%f,%f\n", size, rank, type, runtime, calctime, mpitime, waittime, iotime);
  }

  MPI_Finalize();
  exit(0);
}


void calc(int *iterations, int width, int height, int xstart, int ystart, double xmin, double xmax, double ymin, double ymax, int maxiter) {
  double dx,dy,x,y;
  int    ix,iy;

  dx = (xmax - xmin) / width;
  dy = (ymax - ymin) / height;

  y = ymin;
  for (iy = ystart; iy < ystart + height; ++iy) {
    x = xmin;
    for (ix = xstart; ix < xstart + width; ix++) {
      double zx=0.0,zy=0.0,zxnew;
      int count = 0;
      while ( zx*zx+zy*zy < 16*16 && count < maxiter ) {
        /* z = z*z + (x + i y) */
        zxnew = zx*zx-zy*zy + x;
        zy    = 2*zx*zy     + y;
        zx    = zxnew;
        ++count;
      }
      iterations[iy*width+ix] = count;
      x += dx;
    }
    y += dy;
  }
}
