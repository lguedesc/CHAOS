#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "convergence_test.h"

static void read_params_and_IC(char *name, int *dim, int *npar, int *n, double *t, double *tf, int *ntries, double **par, double **x);

void convergence_test(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *)) {
    
    // Declare Program Parameters
    const double pi = 4 * atan(1);  // Pi number definition
    int DIM;                        // Dimension of the system
    int nPar;                       // Number of parameters of the system
    int nTries;                     // Number of different steps tried
    int N;                          // Number of steps
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;
    double tf;
    double *x = NULL;
    double *par = NULL;
    read_params_and_IC(input_filename, &DIM, &nPar, &N, &t, &tf, &nTries, &par, &x);
    
    /*
    // Time variables
    double time_spent = 0.0;
    clock_t time_i = clock();
    */
    // Call solution
    convergence_test_solution(DIM, N, t, tf, x, nTries, par, edosys);

    /*
    clock_t time_f = clock();
    time_spent += (double)(time_f - time_i) / CLOCKS_PER_SEC; 
    printf("The elapsed time is %f seconds\n", time_spent);
    */
    // Free allocated memory
    free(input_filename);
    free(x); free(par);
}

static void read_params_and_IC(char *name, int *dim, int *npar, int *n, double *t, double *tf, int *ntries, double **par, double **x) {
   // Open input file
    FILE *input = fopen(name, "r");
    if (input == NULL) {
        // Return error if input does not exist 
        perror(name);
        exit(1);
    }
    // Read and assign system constants
    fscanf(input, "%d", dim);
    fscanf(input, "%d", npar);
    // Read and assign program parameters
    fscanf(input, "%d", n); 
    fscanf(input, "%d", ntries);
    // Read and assign initial and final time
    fscanf(input, "%lf %lf", t, tf);
    // Allocate memory for x[dim] and par[npar] vectors
    //*x = malloc((*dim) * sizeof(double));
    *x = malloc((*dim) * sizeof **x);
    *par = malloc((*npar) * sizeof **par);
    // Security check for pointers
    if(*x == NULL || *par == NULL) {
        free(*x); free(*par);
        printf("Memory allocation for *x or *par did not complete successfully");
        return;
    }
    // assign IC to x[dim] vector
    for (int i = 0; i < *dim; i++) {
        fscanf(input, "%lf ", &(*x)[i]);     
    }
    // Assign parameter values to par[npar] vector
    for (int i = 0; i < *npar; i++) {
            fscanf(input, "%lf\n", &(*par)[i]);
    }
    // Close input file
    fclose(input);
    /* The user is responsible to free (x) and (par) after the function call */
}
