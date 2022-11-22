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
#include "epbasin.h"

static void read_params_and_IC(char *name, int *dim, int *npar,  int *np, int *ndiv, double *t, double **par, double **icrange, int *indexX, int *indexY, double **x);
static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, double t, double *x, double *par, double *icrange, int indexX, int indexY, char* funcname, size_t maxlength, double percname, char* mode);

void epbasin(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *)) {
    
    // Parameters related to printing information
    size_t maxLen = 71;             // Max length of the info printed on the screen and on info file
    double percName = 0.6;          // Percentage of space occuped by the name of the quantity printed
    // Declare Program Parameters
    const double pi = 4 * atan(1);  // Pi number definition
    int DIM;                        // Dimension of the system
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int nPar;                       // Number of parameters of the system
    int indexX;                     // Index of control parameter in X
    int indexY;                     // Index of control parameter in Y
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;
    double *x = NULL;
    double *par = NULL;
    double *icRange = NULL;
    read_params_and_IC(input_filename, &DIM, &nPar, &nP, &nDiv, &t, &par, &icRange, &indexX, &indexY, &x);
    
    // Create output files to store results
    char output_epbasin_name[200];
    char output_info_name[200];
    const char *rawdir = "EPBasin/out/";                                                            // Directory of output file
    char *dir = convert_dir(rawdir);
    const char *ext = ".csv";                                                                           // Extension of output file
    const char *ext_info = ".txt";                                                                      // Extension of info file
    snprintf(output_epbasin_name, sizeof(output_epbasin_name), "%s%s_epbasin", dir, outputname);            // Assign name for output file without extension
    snprintf(output_info_name, sizeof(output_info_name), "%s%s_info", dir, outputname);                   // Assign name for output info without extension
    FILE *output_epbasin = create_output_file(output_epbasin_name, ext, dir);                             // Create dynamical diagram output file 
    FILE *output_info = create_output_file(output_info_name, ext_info, dir);                            // Create info output file
    
    // Print information in screen and info output file
    print_info(output_info, DIM, nPar, nP, nDiv, t, x, par, icRange, indexX, indexY, funcname, maxLen, percName, "screen");
    print_info(output_info, DIM, nPar, nP, nDiv, t, x, par, icRange, indexX, indexY, funcname, maxLen, percName, "file");
    
    // Call solution
    ep_basin_of_attraction_2D(output_epbasin, output_info, DIM, nP, nDiv, t, &x, indexX, indexY, icRange, par, nPar, edosys, p_write_epbasin_results);    
    // Close output file
    fclose(output_epbasin);
    fclose(output_info);

    // Free allocated memory
    free(dir);
    free(input_filename);
    free(x); free(par); free(icRange);

}

static void read_params_and_IC(char *name, int *dim, int *npar,  int *np, int *ndiv, double *t, double **par, double **icrange, int *indexX, int *indexY, double **x) {
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
    fscanf(input, "%d %d", np, ndiv); 
    // Read and assign initial time
    fscanf(input, "%lf", t);
    // Allocate memory for x[dim] and par[npar] vectors
    *x = malloc((*dim) * sizeof **x);
    *par = malloc((*npar) * sizeof **par);
    *icrange = malloc (6 * sizeof *icrange);
    // Security check for pointers
    if(*x == NULL || *par == NULL || *icrange == NULL) {
        free(*x); free(*par); free(*icrange);
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
