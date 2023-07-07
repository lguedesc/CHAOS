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
#include "HOS_ftime_series.h"

static void read_params_and_IC(char *name, int *dim, int *npar, int *maxper, int *np, int *ndiv, int* trans, int *nrms, double *t, double **par, double **x, int **rmsindex,
                               int *ncustomvalues, int *nprintf, int *nprintscr, int **printfindex, int **printscrindex);
static void print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, int nrms, double h, double t, double *x, double *par, int *rmsindex, char* funcname, 
                       int ncustomvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex, size_t maxlength, double percname, char* mode);

void HOS_ftime_series(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, size_t, double *, int)) {
    
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
    int attractor;                  // Value of type of motion of the system
    int maxPer;                     // Maximum periodicity to be classified
    int nRMS = 0;                       // Number of state variables that will be submitted to RMS calculation
    int nCustomValues = 0;          // If there is a custom function to be called, this is the number of calculations the function is going to perform
    int nPrintf = 0;                // Number of custom values to be printed in the output file
    int nPrintscr = 0;              // Number of custom values to be printed on the screen
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;
    double *x = NULL;
    double *par = NULL;
    int *rmsindex = NULL;           // Indexes of state variables that will be submitted to RMS calculation
    double *xRMS = NULL;            // State Variables RMS values at permanent regime
    double *overallxRMS = NULL;     // State Variables RMS values at transient + permanent regime
    double *customValues = NULL;    // Variable to store all custom values that custom functions calculate
    char **customNames = NULL;          // Names of the custom values
    int *printfindex = NULL;            // Indexes of custom values that will be printed in the output file
    int *printscrindex = NULL;          // Indexes of custom values that will be printed on the screen
    double *xmin = NULL;                // State variables minimum value at steady state regime
    double *xmax = NULL;                // State variables maximum value at steady state regime
    double *overallxmin = NULL;         // Overall state variables minimum value at steady state regime 
    double *overallxmax = NULL;         // Overall state variables maximum value at steady state regime
    read_params_and_IC(input_filename, &DIM, &nPar, &maxPer, &nP, &nDiv, &trans, &nRMS, &t, &par, &x, &rmsindex, &nCustomValues, &nPrintf, &nPrintscr, &printfindex, &printscrindex);
    // Define Timestep
    double h = (2 * pi) / (nDiv * par[0]); // par[0] = OMEGA
    // Create output files to store results
    char output_ftimeseries_name[200];
    char output_poinc_name[200];
    char output_info_name[200];
    const char *rawdir = "data/FTimeSeries/out/";                                                                // Directory of output file
    char *dir = convert_dir(rawdir);
    const char *ext = ".csv";                                                                               // Extension of output file
    const char *ext_info = ".txt";                                                                          // Extension of info file
    snprintf(output_ftimeseries_name, sizeof(output_ftimeseries_name), "%s%s_ftimeseries", dir, outputname);  // Assign name for output ftimeseries file without extension
    snprintf(output_poinc_name, sizeof(output_poinc_name), "%s%s_poinc", dir, outputname);                    // Assign name for output poincare file without extension
    snprintf(output_info_name, sizeof(output_info_name), "%s%s_info", dir, outputname);                       // Assign name for output info without extension
    FILE *output_ftimeseries = create_output_file(output_ftimeseries_name, ext, dir);                       // Create ftimeseries output file 
    FILE *output_poinc = create_output_file(output_poinc_name, ext, dir);                                   // Create poincare output file 
    FILE *output_info = create_output_file(output_info_name, ext_info, dir);                                     // Create info output file
    
    // Print information in screen and in info file
    print_info(output_info, DIM, nPar, maxPer, nP, nDiv, trans, nRMS, h, t, x, par, rmsindex, funcname, nCustomValues, nPrintf, printfindex, nPrintscr, printscrindex, maxLen, percName, "screen");
    print_info(output_info, DIM, nPar, maxPer, nP, nDiv, trans, nRMS, h, t, x, par, rmsindex, funcname, nCustomValues, nPrintf, printfindex, nPrintscr, printscrindex, maxLen, percName, "file");
    // Call solution
    HOS_full_timeseries_solution(output_ftimeseries, output_poinc, DIM, nP, nDiv, trans, &attractor, maxPer, t, &x, h, par, nRMS, rmsindex, &xRMS, &overallxRMS, &xmin, &xmax, 
                                &overallxmin, &overallxmax, edosys, nCustomValues, &customNames, &customValues, nPrintf, printfindex, nPrintscr, printscrindex, customfunc);
    // Print min and max values on screen and in info file
    print_minmax(xmin, xmax, overallxmin, overallxmax, DIM, maxLen, percName);
    fprint_minmax(output_info, xmin, xmax, overallxmin, overallxmax, DIM, maxLen, percName);
    // Print RMS calculations in screen and in info file
    if (nRMS > 0) {
        print_RMS(nRMS, rmsindex, xRMS, overallxRMS, maxLen, percName);
        fprint_RMS(output_info, nRMS, rmsindex, xRMS, overallxRMS, maxLen, percName);
    }
    // Print custom calculations on screen and in info file
    if (nCustomValues > 0) {
        print_customcalc(nPrintscr, printscrindex, customValues, customNames, maxLen, percName);
        fprint_customcalc(output_info, nPrintscr, printscrindex, customValues, customNames, maxLen, percName);
    }
    // Print attractor in screen and in info file
    print_attractor(attractor, maxPer, maxLen, percName);
    fprint_attractor(output_info ,attractor, maxPer, maxLen, percName);
    // Close output file
    fclose(output_ftimeseries);
    fclose(output_poinc);
    fclose(output_info);

    // Free allocated memory
    free(dir);
    free(input_filename);
    free(x); free(par);
    free(xRMS); free(overallxRMS); free(rmsindex);
    free(xmin); free(xmax); free(overallxmin); free(overallxmax);
    if (nCustomValues > 0) {
        free(customValues); free(printfindex); free(printscrindex);
        for (int i = 0; i < nCustomValues; i++) {
            free(customNames[i]);
        }
        free(customNames); 
    }
}

static void read_params_and_IC(char *name, int *dim, int *npar, int *maxper, int *np, int *ndiv, int* trans, int *nrms, double *t, double **par, double **x, int **rmsindex,
                               int *ncustomvalues, int *nprintf, int *nprintscr, int **printfindex, int **printscrindex) {
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
    fscanf(input, "%d", maxper);
    // Read and assign program parameters
    fscanf(input, "%d %d %d", np, ndiv, trans); 
    // Read and assign initial time
    fscanf(input, "%lf", t);
    // Allocate memory for x[dim] and par[npar] vectors
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
    // Assign number of state variables that will be submitted to RMS calculation
    fscanf(input, "%d", nrms);
    if (nrms > 0) {
        // Allocate memory for rmsindex[nrms]
        *rmsindex = malloc((*nrms) * sizeof **rmsindex);
        // Assign indexes to rmsindex[nrms] vector
        for (int i = 0; i < *nrms; i++) {
                fscanf(input, "%d\n", &(*rmsindex)[i]);
        }
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
        // Assign the number of custom variables that will be printed on screen
        fscanf(input, "%d\n", nprintscr);
        // Allocate memory for printscrindex[nprintf]
        *printscrindex = malloc((*nprintscr) * sizeof **printscrindex);
        // Assign indexes to printscrindex[nprintscr]
        for (int i = 0; i < *nprintscr; i++) {
            fscanf(input, "%d\n", &(*printscrindex)[i]);
        }
    }
    // Close input file
    fclose(input);
    /* The user is responsible to free (x) and (par) after the function call */
}

static void print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, int nrms, double h, double t, double *x, double *par, int *rmsindex, char* funcname, 
                       int ncustomvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex, size_t maxlength, double percname, char* mode) {
    
    if (strcmp(mode, "screen") == 0) {   
        write_prog_parameters_ftimeseries(dim, npar, maxper, np, ndiv, trans, h, maxlength, percname);
        write_initial_conditions(dim, x, t, maxlength, percname);
        write_sys_parameters(npar, par, maxlength, percname);
        write_RMS_calculations_info(nrms, rmsindex, maxlength, percname);
        write_custom_info_calculations(ncustomvalues, nprintf, printfindex, nprintscr, printscrindex, maxlength, percname);
        partition(2, maxlength);
    } 
    else if (strcmp(mode, "file") == 0) {
        fwrite_prog_parameters_ftimeseries(info, funcname, dim, npar, maxper, np, ndiv, trans, h, maxlength, percname);
        fwrite_initial_conditions(info, dim, x, t, maxlength, percname);
        fwrite_sys_parameters(info, npar, par, maxlength, percname);
        fwrite_RMS_calculations_info(info, nrms, rmsindex, maxlength, percname);
        fwrite_custom_info_calculations(info, ncustomvalues, nprintf, printfindex, nprintscr, printscrindex, maxlength, percname);
        fpartition(info, 2, maxlength);
    }
    else {
        printf("Information could not be printed using mode (%s)...\n", mode);
        return;
    }
}
