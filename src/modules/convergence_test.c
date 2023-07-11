#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "../libs/defines.h"
#include "../libs/basic.h"
#include "convergence_test.h"

static void read_params(int dim, int npar, int *n, double *t, double *tf, int *ntries, double **par, double **x);

void convergence_test(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *)) {
    
    // Declare Program Parameters
    int nTries;                     // Number of different steps tried
    int N;                          // Number of steps
    // Assign values for program parameters, system parameters and initial conditions
    double t;
    double tf;
    double *x = NULL;
    double *par = NULL;
    read_params(DIM, nPar, &N, &t, &tf, &nTries, &par, &x);
    // Call solution
    convergence_test_solution(DIM, N, t, tf, x, nTries, par, edosys);
    // Free allocated memory
    free_mem(x, par, NULL);
}

static void read_params(int dim, int npar, int *n, double *t, double *tf, int *ntries, double **par, double **x) {
    // Open input file
    char *input_filename = get_input_filename();
    FILE *input = fopen(input_filename, "r");
    file_safety_check(input);
    // Read and assign program parameters
    fscanf(input, "%d", n); 
    fscanf(input, "%d", ntries);
    // Read and assign initial and final time
    fscanf(input, "%lf %lf", t, tf);
    // Allocate memory for x[dim] and par[npar] vectors
    //*x = malloc((*dim) * sizeof(double));
    *x = malloc(dim * sizeof **x);
    *par = malloc(npar * sizeof **par);
    ptr_safety_check(x, "*x in read_params()");
    ptr_safety_check(par, "*par in read_params()");
    // assign IC to x[dim] vector
    for (int i = 0; i < dim; i++) {
        fscanf(input, "%lf ", &(*x)[i]);     
    }
    // Assign parameter values to par[npar] vector
    for (int i = 0; i < npar; i++) {
            fscanf(input, "%lf\n", &(*par)[i]);
    }
    // Close input file
    fclose(input);
    // Free memory
    free(input_filename);
    /* The user is responsible to free (x) and (par) after the function call */
}
