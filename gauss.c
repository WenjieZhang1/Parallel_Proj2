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

char *ID;

/* Program Parameters */
#define MAXN 6000  /* Max value of N */
int N;  /* Matrix size */
int procs;  /* Number of processors to use */
int myid;
int NTHREADS=4;

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
void gauss1();
/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
  int submit = 0;  /* = 1 if submission parameters should be used */
  int seed = 0;  /* Random seed */
  //char uid[L_cuserid + 2]; /*User name */
  char uid[32];
  /* Read command-line arguments */
  //  if (argc != 3) {
  if ( argc == 1 && !strcmp(argv[1], "submit") ) {
    /* Use submission parameters */
    submit = 1;
    N = 4;
    procs = 2;
    printf("\nSubmission run for \"%s\".\n", uid);
      /*uid = ID;*/
    strcpy(uid,ID);
    srand(randm());
  }
  else {
    if (argc == 3) {
      seed = atoi(argv[3]);
      srand(seed);
      printf("Random seed = %i\n", seed);
    }
    else {
      printf("Usage: %s <matrix_dimension> <num_procs> [random seed]\n",
	     argv[0]);
      printf("       %s submit\n", argv[0]);
      exit(0);
    }
  }
    //  }
  /* Interpret command-line args */
  if (!submit) {
    N = atoi(argv[1]);
    if (N < 1 || N > MAXN) {
      printf("N = %i is out of range.\n", N);
      exit(0);
    }
    procs = atoi(argv[2]);
    if (procs < 1) {
      printf("Warning: Invalid number of processors = %i.  Using 1.\n", procs);
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
  /* Timing variables */
  struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
  struct timezone tzdummy;
  clock_t etstart2, etstop2;  /* Elapsed times using times() */
  unsigned long long usecstart, usecstop;
  struct tms cputstart, cputstop;  /* CPU times for my processes */

  ID = argv[argc-1];
  argc--;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  /* Process program parameters */
  parameters(argc, argv);


  if(myid == 0) {
  	/* Initialize A and B */
    initialize_inputs();
    /* Print input matrices */
    print_inputs();

  /* Start Clock */
  printf("\nStarting clock.\n");
  gettimeofday(&etstart, &tzdummy);
  etstart2 = times(&cputstart);
  }
  /* Gaussian Elimination */
  gauss();
   
  if(myid == 0) {
    /* Stop Clock */
    gettimeofday(&etstop, &tzdummy);
    etstop2 = times(&cputstop);
    printf("Stopped clock.\n");
    usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
    usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

    /* Display output */
    print_X();

    /* Display timing results */
    printf("\nElapsed time = %g ms.\n",
    (float)(usecstop - usecstart)/(float)1000);
  }
  MPI_Finalize();

/* Compare the result
  int right = 1;
  for(int i = 0; i < N; i++) {
	  float dif = X[i] - X1[i];
	  if (dif < 0)	dif = -dif;
	  if(dif > 0.00001) {
  printf("X: %f\n",X[i]);
  printf("X1: %f\n",X1[i]);
		  
		  right = 0;
		  break;
	  }
  }
  printf("right: %d\n",right);
  if(right == 1)	printf("\nRight!\n");
  else	printf("\nWrong!\n");
  */
}


/*----------------Follwoing code used for Testing----------------------------*/
void gauss1() {
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
  MPI_Request request;
  int norm, row, col, i;
  float multiplier;

  for(norm = 0; norm < N - 1; norm++) {
    MPI_Bcast(&A[norm][0], N, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&B[norm], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(myid == 0) {
      for(i = 1; i < procs; i++) {
        for(row = norm+1+i; row < N; row +=procs) {
          MPI_Isend(&A[row], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &request);
          MPI_Wait(&request, &status);
          MPI_Isend(&B[row], 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &request);
          MPI_Wait(&request, &status);
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
          MPI_Recv(&A[row], N, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &status);
          MPI_Recv(&B[row], 1, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &status);
        }
      }
    }
    else {
      for(row = norm + 1 + myid; row < N; row += procs){
        MPI_Recv(&A[row], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&B[row], 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
        /*Gaussian elimination*/
        multiplier = A[row][norm] / A[norm][norm];
        for(col = norm; col < N; col++) {
          A[row][col] -= A[norm][col] * multiplier;
        }
        B[row] -= B[norm] * multiplier;
        MPI_Isend(&A[row], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        MPI_Isend(&B[row], 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if(myid == 0) {
    /* Back substitution */
    for(row = N - 1; row >= 0; row--) {
      X[row] = B[row];
      for(col = N-1; col > row; col--) {
        X[row] -= A[row][col] * X[col];
      }
      X[row] /= A[row][row];
    }
  }
}
