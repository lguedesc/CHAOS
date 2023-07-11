#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "../libs/interface.h"
#include "../libs/defines.h"
#include "../libs/basic.h"
#include "epbasin.h"

static void read_params(int dim, int npar,  int *np, int *ndiv, double *t, double **par, double **icrange, int *indexX, int *indexY, double **x);
static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, double t, double *x, double *par, double *icrange, int indexX, int indexY, char* funcname, size_t maxlength, double percname, char* mode);

void epbasin(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *)) {
    
    // Declare Program Parameters
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int indexX;                     // Index of control parameter in X
    int indexY;                     // Index of control parameter in Y
    // Assign values for program parameters, system parameters and initial conditions
    double t;
    double *x = NULL;
    double *par = NULL;
    double *icRange = NULL;
    read_params(DIM, nPar, &nP, &nDiv, &t, &par, &icRange, &indexX, &indexY, &x);
    // Create output files to store results
    const char *directory = "data/EPBasin/out/";                                                    // Directory of output file
    const char *module = "epbasin";
    FILE *output_epbasin = name_and_create_output_files(outputname, directory, module, ".csv");     // Create dynamical diagram output file 
    FILE *output_info = name_and_create_output_files(outputname, directory, "info", ".txt");        // Create info output file
    // Print information in screen and info output file
    print_info(output_info, DIM, nPar, nP, nDiv, t, x, par, icRange, indexX, indexY, funcname, MAX_PRINT_LEN, PERC_PRINT_NAME, "screen");
    print_info(output_info, DIM, nPar, nP, nDiv, t, x, par, icRange, indexX, indexY, funcname, MAX_PRINT_LEN, PERC_PRINT_NAME, "file");
    // Call solution
    ep_basin_of_attraction_2D(output_epbasin, output_info, DIM, nP, nDiv, t, &x, indexX, indexY, icRange, par, nPar, edosys, p_write_epbasin_results);    
    // Close output file
    close_files(2, output_epbasin, output_info);
    // Free allocated memory
    free_mem(x, par, icRange, NULL);
}

static void read_params(int dim, int npar,  int *np, int *ndiv, double *t, double **par, double **icrange, int *indexX, int *indexY, double **x) {
    // Open input file
    char *input_filename = get_input_filename();
    FILE *input = fopen(input_filename, "r");
    file_safety_check(input);
    // Read and assign program parameters
    fscanf(input, "%d %d", np, ndiv); 
    // Read and assign initial time
    fscanf(input, "%lf", t);
    // Allocate memory for x[dim] and par[npar] vectors
    *x = malloc(dim * sizeof **x);
    *par = malloc(npar * sizeof **par);
    *icrange = malloc (6 * sizeof *icrange);
    ptr_safety_check(x, "*x in read_params()");
    ptr_safety_check(par, "*par in read_params()");
    ptr_safety_check(icrange, "*icrange in read_params()");
    // assign IC to x[dim] vector
    for (int i = 0; i < dim; i++) {
        fscanf(input, "%lf ", &(*x)[i]);     
    }
    // Assign parameter values to par[npar] vector
    for (int i = 0; i < npar; i++) {
        fscanf(input, "%lf\n", &(*par)[i]);
    }
    // Assign index of X control parameter of the diagram
    fscanf(input, "%d\n", indexX);
    // Assign range of control parameter in X
    for (int i = 0; i < 3; i++) {
        fscanf(input, "%lf\n", &(*icrange)[i]);
    }
    // Assign index of Y control parameter of the diagram
    fscanf(input, "%d\n", indexY);
    // Assign range of control parameter in Y
    for (int i = 3; i < 6; i++) {
        fscanf(input, "%lf\n", &(*icrange)[i]);
    }
    // Close input file
    fclose(input);
    // Free memory
    free(input_filename);
    /* The user is responsible to free (x), (par) or (parrange) after the function call */
}

static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, double t, double *x, double *par, double *icrange, int indexX, int indexY, char* funcname, size_t maxlength, double percname, char* mode) {

    char *type = "Fixed Points";
    if (strcmp(mode, "screen") == 0) {
        write_prog_parameters_epbasin(dim, npar, np, ndiv, maxlength, percname);
        write_initial_conditions(dim, x, t, maxlength, percname);
        write_sys_parameters(npar, par, maxlength, percname);
        write_basin_info(type, icrange, indexX, indexY, maxlength, percname);
        partition(2, maxlength);          
    } 
    else if (strcmp(mode, "file") == 0) {
        fwrite_prog_parameters_epbasin(info, funcname, dim, npar, np, ndiv, maxlength, percname);
        fwrite_initial_conditions(info, dim, x, t, maxlength, percname);
        fwrite_sys_parameters(info, npar, par, maxlength, percname);
        fwrite_basin_info(type, info, icrange, indexX, indexY, maxlength, percname);
        fpartition(info, 2, maxlength);
    }
    else {
        printf("  Information could not be printed in file using mode (%s)...\n", mode);
        return;
    }

}
