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
// void gauss_test();
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
  if(argc != 1) {
    printf("\n Please only enter the matrix size");
  }
  else {
    seed = 5555;
    srand(seed);
    N = atoi(argv[1]);
    if(N<1 || N>MAXN) {
      printf("N = %i is out of range.\n", N);
      exit(0);
    }
    /* Print parameters */
    if(id == 0) {
      printf("\nMatrix dimension N = %i.\n", N);
    }
  }
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
  argc--;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  /* Process program parameters */
  parameters(argc, argv);

  /* Run Gaussian Function */
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
    if(right == 1)  printf("\nThe Result is Right!\n");
    else  printf("\nThe Result is Wrong!\n");
  }

  MPI_Finalize();
  return 0;
}


/*----------------Follwoing code used for Testing----------------------------*/
void gauss_test() {
  int norm, row, col; 
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


void gauss() {
  int norm, row, col;
  float multiplier;
  double startwtime = 0.0;
  double endwtime;

  if(id == 0) {
    initialize_inputs();
    /* Print input matrices */
    print_inputs();
    printf("\nStart Computing Parallely Using MPI.\n");
    startwtime = MPI_Wtime();
  }
  MPI_Bcast(&A, N*N, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&B, N, MPI_FLOAT, 0, MPI_COMM_WORLD);
  for(norm = 0; norm < N - 1; norm++) {
      /*Gaussian elimination*/
      row = id+1;
      while(row < norm){
        row += procs;
      }
      for(; row < N; row += procs) {
        multiplier = A[row][norm] / A[norm][norm];
        for(col = norm; col < N; col++) {
          A[row][col] -= A[norm][col] * multiplier;
        }
        B[row] -= B[norm] * multiplier;
      }
      MPI_Bcast(&A[norm+1], N, MPI_FLOAT, norm%procs, MPI_COMM_WORLD);
      MPI_Bcast(&B[norm+1], 1, MPI_FLOAT, norm%procs, MPI_COMM_WORLD);
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
