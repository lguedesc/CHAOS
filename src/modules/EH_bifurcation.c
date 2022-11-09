#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "../libs/energyharvest.h"
#include "../libs/interface.h"
#include "EH_bifurcation.h"

static void read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *parindex, double **x, int *nrms, int **rmsindex, int *bifmode,
                                int *ncustomvalues, int *nprintf, int **printfindex);
static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int parindex, int nrms, int *rmsindex, int bifmode, char* funcname,
                        int ncustomvalues, int nprintf, int *printfindex, size_t maxlength, double percname, char* mode);

void EH_bifurcation(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, int, int, char **, size_t, double *, int)) {
    // Parameters related to printing information
    size_t maxLen = 71;             // Max length of the info printed on the screen and on info file
    double percName = 0.6;          // Percentage of space occuped by the name of the quantity printed
    // Declare Program Parameters
    const double pi = 4 * atan(1);  // Pi number definition
    int DIM;                        // Dimension of the system
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int nPar;                       // Number of parameters of the system
    int trans;                      // Value of nP in which greater values are considered transient response
    int bMode;                      // Mode of the Bifurcation Diagram: 0 to follow attractor, 1 to reset initial conditions in each step
    int parIndex;                   // Index of control parameter in par array
    int nRMS = 0;                       // Number of state variables that will be submitted to RMS calculation
    int nCustomValues = 0;          // If there is a custom function to be called, this is the number of calculations the function is going to perform
    int nPrintf = 0;                // Number of custom values to be printed in the output file
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;
    double *x = NULL;
    double *par = NULL;
    double *parRange = NULL;
    int *rmsindex = NULL;           // Indexes of state variables that will be submitted to RMS calculation
    int *printfindex = NULL;        // Indexes of custom values that will be printed in the output file
    read_params_and_IC(input_filename, &DIM, &nPar, &nP, &nDiv, &trans, &t, &par, &parRange, &parIndex, &x, &nRMS, &rmsindex, &bMode,
                       &nCustomValues, &nPrintf, &printfindex);
    
    // Create output files to store results
    char output_bifurc_name[200];
    char output_info_name[200];
    char output_poinc_name[200];
    const char *rawdir = "Bifurcation/out/";                                                            // Directory of output file
    char *dir = convert_dir(rawdir);
    const char *ext = ".csv";                                                                           // Extension of output file
    const char *ext_info = ".txt";                                                                      // Extension of info file
    snprintf(output_poinc_name, sizeof(output_poinc_name), "%s%s_bifurc_poinc", dir, outputname);     // Assign name for output rk4 without extension
    snprintf(output_info_name, sizeof(output_info_name), "%s%s_info", dir, outputname);                 // Assign name for output info without extension
    snprintf(output_bifurc_name, sizeof(output_bifurc_name), "%s%s_bifurc", dir, outputname);           // Assign name for output rk4 without extension
    FILE *output_bifurc_poinc = create_output_file(output_poinc_name, ext, dir);                        // Create poincare bifurc output file 
    FILE *output_bifurc = create_output_file(output_bifurc_name, ext, dir);                             // Create bifurc output file
    FILE *output_info = create_output_file(output_info_name, ext_info, dir);                            // Create info output file
    
    // Print information in screen and info output file
    print_info(output_info, DIM, nPar, nP, nDiv, trans, t, x, par, parRange, parIndex, nRMS, rmsindex, bMode, funcname, nCustomValues, nPrintf, printfindex, maxLen, percName, "screen");
    print_info(output_info, DIM, nPar, nP, nDiv, trans, t, x, par, parRange, parIndex, nRMS, rmsindex, bMode, funcname, nCustomValues, nPrintf, printfindex, maxLen, percName, "file");
    // To store the runtime of the program
    //double time_spent = 0.0;
    //clock_t time_i = clock();
    // Call solution
    EH_bifurc_solution(output_bifurc, output_bifurc_poinc, DIM, nP, nDiv, trans, t, x, parIndex, parRange, par, nRMS, rmsindex, edosys, nCustomValues, nPrintf, printfindex, customfunc, bMode);
    // Close output file
    fclose(output_bifurc);
    fclose(output_info);

    // Free allocated memory
    free(dir);
    free(input_filename);
    free(x); free(par); free(parRange);
    free(rmsindex);
    free(printfindex);

    // Calculate time of execution
    //clock_t time_f = clock();
    //time_spent += (double)(time_f - time_i) / CLOCKS_PER_SEC; 
    //printf("The elapsed time is %f seconds", time_spent);

}

static void read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *parindex, double **x, int *nrms, int **rmsindex, int *bifmode,
                                int *ncustomvalues, int *nprintf, int **printfindex) {
    // Open input file
    FILE *input = fopen(name, "r");
    if (input == NULL) {
        // Return error if input does not exist 
        perror(name);
        exit(1);
    }
    // Determine the mode of the bifurcation diagram: 0 to follow attractor, 1 to reset initial conditions in each step
    fscanf(input, "%d", &(*bifmode));
    // Read and assign system constants
    fscanf(input, "%d", dim);
    fscanf(input, "%d", npar);
    // Read and assign program parameters
    fscanf(input, "%d %d %d", np, ndiv, trans); 
    // Read and assign initial time
    fscanf(input, "%lf", t);
    // Allocate memory for x[dim] and par[npar] vectors
    *x = malloc((*dim) * sizeof **x);
    *par = malloc((*npar) * sizeof **par);
    *parrange = malloc (3 * sizeof *parrange);
    // Security check for pointers
    if(*x == NULL || *par == NULL || *parrange == NULL) {
        free(*x); free(*par); free(*parrange);
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
    // Assign index of control parameter of the diagram
    fscanf(input, "%d\n", parindex);
    // Assign range of parameter of bifurcation
    for (int i = 0; i < 3; i++) {
        fscanf(input, "%lf\n", &(*parrange)[i]);
    }
    // Assign number of state variables that will be submitted to RMS calculation
    fscanf(input, "%d", nrms);
    // Allocate memory for rmsindex[nrms]
    *rmsindex = malloc((*nrms) * sizeof **rmsindex);
    // Assign indexes to rmsindex[nrms] vector
    for (int i = 0; i < *nrms; i++) {
            fscanf(input, "%d\n", &(*rmsindex)[i]);
    }
    // Assign the number of custom variables to be calculated
    fscanf(input, "%d\n", ncustomvalues);
    if (ncustomvalues > 0) {
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
    /* The user is responsible to free (x), (par) or (parrange) after the function call */
}

static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int parindex, int nrms, int *rmsindex, int bifmode, char* funcname,
                        int ncustomvalues, int nprintf, int *printfindex, size_t maxlength, double percname, char* mode) {

    if (strcmp(mode, "screen") == 0) {   
        write_prog_parameters_bifurcation(dim, npar, np, ndiv, trans, maxlength, percname);
        write_initial_conditions(dim, x, t, maxlength, percname);
        write_sys_parameters(npar, par, maxlength, percname);
        write_bifurcation_info(parrange, parindex, bifmode, maxlength, percname);
        write_RMS_calculations_info(nrms, rmsindex, maxlength, percname);
        write_custom_info_calculations(ncustomvalues, nprintf, printfindex, 0, NULL, maxlength, percname);
        partition(2, maxlength);
    } 
    else if (strcmp(mode, "file") == 0) {
        fwrite_prog_parameters_bifurcation(info, funcname, dim, npar, np, ndiv, trans, maxlength, percname);
        fwrite_initial_conditions(info, dim, x, t, maxlength, percname);
        fwrite_sys_parameters(info, npar, par, maxlength, percname);
        fwrite_bifurcation_info(info, parrange, parindex, bifmode, maxlength, percname);
        fwrite_RMS_calculations_info(info, nrms, rmsindex, maxlength, percname);
        fwrite_custom_info_calculations(info, ncustomvalues, nprintf, printfindex, 0, NULL, maxlength, percname);
        fpartition(info, 2, maxlength);
    }
    else {
        printf("Information could not be printed using mode (%s)...\n", mode);
        return;
    }

}
