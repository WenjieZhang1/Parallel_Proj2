/* Gaussian elimination without pivoting.
 */

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
 * You need not submit the provided code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include "mpi.h"


/*#include <ulocks.h>
#include <task.h>
*/


/* Program Parameters */
#define MAXN 6000  /* Max value of N */
int N;  /* Matrix size */
int procs;  /* Number of processors to use */
int id;

/* Matrices and vectors */
volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
volatile float A1[MAXN][MAXN], B1[MAXN], X1[MAXN];
/* A * X = B, solve for X */

/* junk */
#define randm() 4|2[uid]&3

/* Prototype */
void gauss();  /* The function you will provide.
		* It is this routine that is timed.
		* It is called only on the parent.
		*/
void gauss_test();
/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
  int seed = 0;  /* Random seed */

  /* Read command-line arguments */
  if(argc != 3) {
    printf("\n Please type three arguments");
  }
  else {
    seed = atoi(argv[3]);
    srand(seed);
    printf("Random seed =  %i\n", seed);
    N = atoi(argv[1]);
    if(N<1 || N>MAXN) {
      printf("N = %i is out of range.\n", N);
      exit(0);
    }
    procs = atoi(argv[2]);
    if(procs < 1) {
      printf("Warning: Invalid number of processor = %i. Using 1.\n", procs);
      procs = 1;
    }
  }
  /* Print parameters */
  printf("\nMatrix dimension N = %i.\n", N);
  printf("Number of processors = %i.\n", procs);
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
  int row, col;

  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
	  A1[row][col] = A[row][col];
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
	B1[col] = B[col];
	X1[col] = X[col];
  }

}

/* Print input matrices */
void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
	printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
  if (N < 10) {
    printf("\nA1 =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
	printf("%5.2f%s", A1[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB1 = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B1[col], (col < N-1) ? "; " : "]\n");
    }
  }
}

void print_X() {
  int row;

  if (N < 10) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
  if (N < 10) {
    printf("\nX1 = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X1[row], (row < N-1) ? "; " : "]\n");
    }
  }
}

int main(int argc, char **argv) {
  // /* Timing variables */
  // struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
  // struct timezone tzdummy;
  // clock_t etstart2, etstop2;  /* Elapsed times using times() */
  // unsigned long long usecstart, usecstop;
  // struct tms cputstart, cputstop;  /* CPU times for my processes */

  argc--;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  /* Process program parameters */
  parameters(argc, argv);


  if(id == 0) {
  	/* Initialize A and B */
    initialize_inputs();
    /* Print input matrices */
    print_inputs();
  }

  /* Gaussian Elimination */
  gauss();
   
  if(id == 0) {
    /* Display output */
    print_X();
    gauss_test();
    /* Compare the result*/

    int right = 1;
    int j = 0;
    for(; j < N; j++) {
      float dif = X[j] - X1[j];
      if (dif < 0)  dif = -dif;
      if (dif > 0.0001) {
        printf("X: %f\n", X[j]);
        printf("X1: %f\n", X1[j]);
        right = 0;
        break;
      }
    }

    printf("right: %d\n",right);
    if(right == 1)  printf("\nRight!\n");
    else  printf("\nWrong!\n");
  
  }

  MPI_Finalize();
  return 0;
}


/*----------------Follwoing code used for Testing----------------------------*/
void gauss_test() {
  int norm, row, col;  /* Normalization row, and zeroing
			* element row and col */
  float multiplier;

  printf("Computing Serially.\n");

  /* Gaussian elimination */
  for (norm = 0; norm < N - 1; norm++) {
    for (row = norm + 1; row < N; row++) {
      multiplier = A1[row][norm] / A1[norm][norm];
      for (col = norm; col < N; col++) {
		  A1[row][col] -= A1[norm][col] * multiplier;
      }
      B1[row] -= B1[norm] * multiplier;
    }
  }
  /* Back substitution */
  for (row = N - 1; row >= 0; row--) {
    X1[row] = B1[row];
    for (col = N-1; col > row; col--) {
      X1[row] -= A1[row][col] * X1[col];
    }
    X1[row] /= A1[row][row];
  }
}
/*-----------Above code used for testing------------------------*/


/* ------------------ Above Was Provided --------------------- */

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, procs, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */

void gauss() {
  MPI_Status status;
  int norm, row, col, i;
  float multiplier;
  double startwtime = 0.0;
  double endwtime;
  MPI_Barrier(MPI_COMM_WORLD);

  if(id == 0) {
    printf("\nStart Computing Parallely Using MPI.\n");
    startwtime = MPI_Wtime();
  }
  for(norm = 0; norm < N - 1; norm++) {
    MPI_Bcast(&(A[norm][0]), N, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(B[norm]), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(id == 0) {
      for(i = 1; i < procs; i++) {
        for(row = norm+1+i; row < N; row +=procs) {
          MPI_Send(&(A[row]), N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
          MPI_Send(&(B[row]), 1, MPI_FLOAT, i, 1, MPI_COMM_WORLD);
        }
      }
      /*Gaussian elimination*/
      for(row = norm + 1; row < N; row += procs) {
        multiplier = A[row][norm] / A[norm][norm];
        for(col = norm; col < N; col++) {
          A[row][col] -= A[norm][col] * multiplier;
        }
        B[row] -= B[norm] * multiplier;
      }
      /*Receive the updated data from other processes*/
      for(i = 1; i < procs; i++){
        for(row = norm + 1 + i; row < N; row += procs) {
          MPI_Recv(&(A[row]), N, MPI_FLOAT, i, 2, MPI_COMM_WORLD, &status);
          MPI_Recv(&(B[row]), 1, MPI_FLOAT, i, 3, MPI_COMM_WORLD, &status);
        }
      }
    }
    else {
      for(row = norm + 1 + id; row < N; row += procs){
        MPI_Recv(&A[row], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&B[row], 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
        /*Gaussian elimination*/
        multiplier = A[row][norm] / A[norm][norm];
        for(col = norm; col < N; col++) {
          A[row][col] -= A[norm][col] * multiplier;
        }
        B[row] -= B[norm] * multiplier;
        MPI_Send(&A[row], N, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&B[row], 1, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if(id == 0) {
    /* Back substitution */
    for(row = N - 1; row >= 0; row--) {
      X[row] = B[row];
      for(col = N-1; col > row; col--) {
        X[row] -= A[row][col] * X[col];
      }
      X[row] /= A[row][row];
    }
    endwtime = MPI_Wtime();
    printf("elapsed time = %f\n", endwtime - startwtime);
  }
}
