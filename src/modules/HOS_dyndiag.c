#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "../libs/nlosc.h"
#include "../libs/interface.h"
#include "../libs/defines.h"
#include "../libs/basic.h"
#include "HOS_dyndiag.h"

static void read_params(int dim, int npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *indexX, int *indexY, double **x, int *nrms, int **rmsindex, int *bifmode,
                        int *ncustomvalues, int *nprintf, int **printfindex);
static void print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int indexX, int indexY, int nrms, int *rmsindex, int bifmode,int ncustomvalues, int nprintf, int *printfindex, size_t maxlength, double percname, char* funcname, char* mode);

void HOS_dyndiag(char *funcname, unsigned int DIM, unsigned int nPar, ang_info *angles, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int)) {
    // Declare Program Parameters
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int trans;                      // Value of nP in which greater values are considered transient response
    int dMode;                      // Mode of the Dynamical Diagram: 0 to follow attractor, 1 to reset initial conditions in each step
    int maxPer;                     // Maximum periodicity to be classified
    int indexX;                     // Index of control parameter in X
    int indexY;                     // Index of control parameter in Y
    int nRMS = 0;                   // Number of state variables that will be submitted to RMS calculation
    int nCustomValues = 0;          // If there is a custom function to be called, this is the number of calculations the function is going to perform
    int nPrintf = 0;                // Number of custom values to be printed in the output file
    // Assign values for program parameters, system parameters and initial conditions
    double t;
    double *x = NULL;
    double *par = NULL;
    double *parRange = NULL;
    int *rmsindex = NULL;           // Indexes of state variables that will be submitted to RMS calculation
    int *printfindex = NULL;        // Indexes of custom values that will be printed in the output file
    read_params(DIM, nPar, &maxPer, &nP, &nDiv, &trans, &t, &par, &parRange, &indexX, &indexY, &x, &nRMS, &rmsindex, &dMode, 
                &nCustomValues, &nPrintf, &printfindex);
    // Create output files to store results
    const char *directory = "data/DynDiagram/out/";                                                     // Directory of output file
    const char *module = "dyndiag";
    FILE *output_dyndiag = name_and_create_output_files(outputname, directory, module, ".csv");         // Create dynamical diagram output file 
    FILE *output_info = name_and_create_output_files(outputname, directory, "info", ".txt");            // Create info output file
    // Print information in screen and info output file
    print_info(output_info, DIM, nPar, maxPer, nP, nDiv, trans, t, x, par, parRange, indexX, indexY, nRMS, rmsindex, dMode, nCustomValues, nPrintf, printfindex, MAX_PRINT_LEN, PERC_PRINT_NAME, funcname, "screen");
    print_info(output_info, DIM, nPar, maxPer, nP, nDiv, trans, t, x, par, parRange, indexX, indexY, nRMS, rmsindex, dMode, nCustomValues, nPrintf, printfindex, MAX_PRINT_LEN, PERC_PRINT_NAME, funcname, "file");
    // Call solution
    HOS_dynamical_diagram_solution(output_dyndiag, DIM, nP, nDiv, trans, maxPer, t, &x, indexX, indexY, parRange, par, angles, nPar, nRMS, rmsindex, edosys,
                                    nCustomValues, nPrintf, printfindex, customfunc, dMode);    
    // Close output file
    close_files(2, output_dyndiag, output_info);
    // Free allocated memory
    free_mem(x, par, parRange, rmsindex, printfindex, NULL);
}

static void read_params(int dim, int npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *indexX, int *indexY, double **x, int *nrms, int **rmsindex, int *bifmode,
                        int *ncustomvalues, int *nprintf, int **printfindex) {
    // Open input file
    char *input_filename = get_input_filename();
    FILE *input = fopen(input_filename, "r");
    file_safety_check(input);
    // Determine the mode of the bifurcation diagram: 0 to follow attractor, 1 to reset initial conditions in each step
    fscanf(input, "%d", &(*bifmode));
    // Read and assign system constants
    fscanf(input, "%d", maxper);
    fscanf(input, "%d %d %d", np, ndiv, trans); 
    // Read and assign initial time
    fscanf(input, "%lf", t);
    // Allocate memory for x[dim] and par[npar] vectors
    *x = malloc(dim * sizeof **x);
    *par = malloc(npar * sizeof **par);
    *parrange = malloc (6 * sizeof *parrange);
    ptr_safety_check(x, "*x in read_params()");
    ptr_safety_check(par, "*par in read_params()");
    ptr_safety_check(parrange, "*parrange in read_params()");
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
        fscanf(input, "%lf\n", &(*parrange)[i]);
    }
    // Assign index of Y control parameter of the diagram
    fscanf(input, "%d\n", indexY);
    // Assign range of control parameter in Y
    for (int i = 3; i < 6; i++) {
        fscanf(input, "%lf\n", &(*parrange)[i]);
    }
    // Assign number of state variables that will be submitted to RMS calculation
    fscanf(input, "%d", nrms);
    // Allocate memory for rmsindex[nrms]
    *rmsindex = malloc((*nrms) * sizeof **rmsindex);
    // Assign indexes to rmsindex[nrms] vector
    for (int i = 0; i < (*nrms); i++) {
            fscanf(input, "%d\n", &(*rmsindex)[i]);
    }
    // Assign the number of custom variables to be calculated
    fscanf(input, "%d\n", ncustomvalues);
    if ((*ncustomvalues) > 0) {
        // Assign the number of custom variables that will be printed in output file
        fscanf(input, "%d\n", nprintf);
        // Allocate memory for printfindex[nprintf]
        *printfindex = malloc((*nprintf) * sizeof **printfindex);
        // Assign indexes to printfindex[nprintf]
        for (int i = 0; i < *nprintf; i++) {
            fscanf(input, "%d\n", &(*printfindex)[i]);
        }
    }
    // Close input file
    fclose(input);
    // Free memory
    free(input_filename);
    /* The user is responsible to free everything allocated after the function call */
}

static void print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int indexX, int indexY, int nrms, int *rmsindex, int bifmode, int ncustomvalues, int nprintf, int *printfindex, size_t maxlength, double percname, char* funcname, char* mode) {
    if (strcmp(mode, "screen") == 0) {   
        write_prog_parameters_dyndiag(dim, npar, np, ndiv, trans, maxlength, percname);
        write_initial_conditions(dim, x, t, maxlength, percname);
        write_sys_parameters(npar, par, maxlength, percname);
        write_dyndiag_info(parrange, indexX, indexY, bifmode, maxlength, percname);
        write_RMS_calculations_info(nrms, rmsindex, maxlength, percname);
        write_custom_info_calculations(ncustomvalues, nprintf, printfindex, 0, NULL, maxlength, percname);
        partition(2, maxlength);        
    } 
    else if (strcmp(mode, "file") == 0) {
        fwrite_prog_parameters_dyndiag(info, funcname, dim, npar, np, ndiv, trans, maxlength, percname);
        fwrite_initial_conditions(info, dim, x, t, maxlength, percname);
        fwrite_sys_parameters(info, npar, par, maxlength, percname);
        fwrite_dyndiag_info(info, parrange, indexX, indexY, bifmode, maxlength, percname);
        fwrite_RMS_calculations_info(info, nrms, rmsindex, maxlength, percname);
        fwrite_custom_info_calculations(info, ncustomvalues, nprintf, printfindex, 0, NULL, maxlength, percname);
        fpartition(info, 2, maxlength);
    }
    else {
        printf("Information could not be printed in file using mode (%s)...\n", mode);
        return;
    }

}
