#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../libs/interface.h"
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "../libs/nlosc.h" 
#include "../libs/defines.h"
#include "../libs/basic.h"
#include "HOS_timeseries.h"

static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, int nrms, double h, double t, double *x, double *par, int *rmsindex, char* funcname, 
                       int ncustomvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex, size_t maxlength, double percname, char* mode);
static void read_params(int dim, int npar, int *np, int *ndiv, int *trans, int *nrms, double *t, double **par, double **x, int **rmsindex, int *ncustomvalues, int *nprintf, int *nprintscr, int **printfindex, int **printscrindex);

static void print_results(FILE *info, int DIM, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int nRMS, int *rmsindex, double *xRMS, double *overallxRMS,
                          int nCustomValues, int nPrintscr, int *printscrindex, double *customValues, char **customNames);

void HOS_timeseries(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int)) {
    // Declare Program Parameters 
    int nP;                                     // Number of forcing periods analyzed
    int nDiv;                                   // Number of divisions in each forcing period
    int trans;                                  // Number of forcing periods considered transient
    int nRMS = 0;                               // Number of state variables that will be submitted to RMS calculation
    int nCustomValues = 0;                      // If there is a custom function to be called, this is the number of calculations the function is going to perform
    int nPrintf = 0;                            // Number of custom values to be printed in the output file
    int nPrintscr = 0;                          // Number of custom values to be printed on the screen
    double t;                                   // Time
    double *x = NULL;                           // State Variables
    double *par = NULL;                         // Parameters
    int *rmsindex = NULL;                       // Indexes of state variables that will be submitted to RMS calculation
    double *xRMS = NULL;                        // State Variables RMS values at permanent regime
    double *overallxRMS = NULL;                 // State Variables RMS values at transient + permanent regime
    double *customValues = NULL;                // Variable to store all custom values that custom functions calculate
    char **customNames = NULL;                  // Names of the custom values
    int *printfindex = NULL;                    // Indexes of custom values that will be printed in the output file
    int *printscrindex = NULL;                  // Indexes of custom values that will be printed on the screen
    double *xmin = NULL;                        // State variables minimum value at steady state regime
    double *xmax = NULL;                        // State variables maximum value at steady state regime
    double *overallxmin = NULL;                 // Overall state variables minimum value at steady state regime 
    double *overallxmax = NULL;                 // Overall state variables maximum value at steady state regime 

    read_params(DIM, nPar, &nP, &nDiv, &trans, &nRMS, &t, &par, &x, &rmsindex, &nCustomValues, &nPrintf, &nPrintscr, &printfindex, &printscrindex);
    // Define Timestep
    double h = (2 * PI) / (nDiv * par[0]); // par[0] = OMEGA
    // Create output files
    const char *directory = "data/TimeSeries/out/";
    const char *module = "timeseries";
    FILE *output_timeseries = name_and_create_output_files(outputname, directory, module, ".csv");
    FILE *output_info = name_and_create_output_files(outputname, directory, "info", ".txt");
    // Print information in screen and info output file
    print_info(output_info, DIM, nPar, nP, nDiv, trans, nRMS, h, t, x, par, rmsindex, funcname, nCustomValues, nPrintf, printfindex, nPrintscr, printscrindex, MAX_PRINT_LEN, PERC_PRINT_NAME, "screen");
    print_info(output_info, DIM, nPar, nP, nDiv, trans, nRMS, h, t, x, par, rmsindex, funcname, nCustomValues, nPrintf, printfindex, nPrintscr, printscrindex, MAX_PRINT_LEN, PERC_PRINT_NAME,"file");
    /*
    // Time variables
    double time_spent = 0.0;
    clock_t time_i = clock();
    */
    // Call solution
    HOS_timeseries_solution(output_timeseries, DIM, nP, nDiv, trans, t, x, h, par, nRMS, rmsindex, &xRMS, &overallxRMS, &xmin, &xmax, &overallxmin, &overallxmax, edosys,  
                            nCustomValues, &customNames, &customValues, nPrintf, printfindex, nPrintscr, printscrindex, customfunc);    
    /*
    clock_t time_f = clock();
    time_spent += (double)(time_f - time_i) / CLOCKS_PER_SEC; 
    printf("The elapsed time is %f seconds\n", time_spent);
    */
    
    // Print some results on the screen
    print_results(output_info, DIM, xmin, xmax, overallxmin, overallxmax, nRMS, rmsindex, xRMS, overallxRMS, nCustomValues, nPrintscr, printscrindex, customValues, customNames);   
    // Close output file
    close_files(2, output_timeseries, output_info);
    // Free allocated memory
    free_mem(x, par, xRMS, overallxRMS, rmsindex, xmin, xmax, overallxmin, overallxmax, NULL);
    if (nCustomValues > 0) {
        free_mem(customValues, printfindex, printscrindex, NULL);
        free_2D_mem((void **)customNames, nCustomValues);         
    }
}

static void read_params(int dim, int npar, int *np, int *ndiv, int *trans, int *nrms, double *t, double **par, double **x, int **rmsindex, int *ncustomvalues, int *nprintf, int *nprintscr, int **printfindex, int **printscrindex) {
    // Open input file
    char *input_filename = get_input_filename();
    FILE *input = fopen(input_filename, "r");
    file_safety_check(input);
    // Read and assign program parameters
    fscanf(input, "%d %d %d", np, ndiv, trans); 
    // Read and assign initial time
    fscanf(input, "%lf", t);
    // Allocate memory for x[dim] and par[npar] vectors
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
    // Free memory
    free(input_filename);
    /* The user is responsible to free (x) and (par) after the function call */
}

static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, int nrms, double h, double t, double *x, double *par, int *rmsindex, char* funcname, 
                       int ncustomvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex, size_t maxlength, double percname, char* mode) {
    
    if (strcmp(mode, "screen") == 0) {   
        write_prog_parameters_timeseries(dim, npar, np, ndiv, trans, h, maxlength, percname);
        write_initial_conditions(dim, x, t, maxlength, percname);
        write_sys_parameters(npar, par, maxlength, percname);
        write_RMS_calculations_info(nrms, rmsindex, maxlength, percname);
        write_custom_info_calculations(ncustomvalues, nprintf, printfindex, nprintscr, printscrindex, maxlength, percname);
    } 
    else if (strcmp(mode, "file") == 0) {
        fwrite_prog_parameters_timeseries(info, funcname, dim, npar, np, ndiv, trans, h, maxlength, percname);
        fwrite_initial_conditions(info, dim, x, t, maxlength, percname);
        fwrite_sys_parameters(info, npar, par, maxlength, percname);
        fwrite_RMS_calculations_info(info, nrms, rmsindex, maxlength, percname);
        fwrite_custom_info_calculations(info, ncustomvalues, nprintf, printfindex, nprintscr, printscrindex, maxlength, percname);
    }
    else {
        printf("Information could not be printed using mode (%s)...\n", mode);
        return;
    }
}

static void print_results(FILE *info, int DIM, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int nRMS, int *rmsindex, double *xRMS, double *overallxRMS,
                          int nCustomValues, int nPrintscr, int *printscrindex, double *customValues, char **customNames) {
    // Print min and max values on screen and in info file
    print_minmax(xmin, xmax, overallxmin, overallxmax, DIM, MAX_PRINT_LEN, PERC_PRINT_NAME);
    fprint_minmax(info, xmin, xmax, overallxmin, overallxmax, DIM, MAX_PRINT_LEN, PERC_PRINT_NAME);
    // Print RMS calculations in screen and in info file
    if (nRMS > 0) {
        print_RMS(nRMS, rmsindex, xRMS, overallxRMS, MAX_PRINT_LEN, PERC_PRINT_NAME);
        fprint_RMS(info, nRMS, rmsindex, xRMS, overallxRMS, MAX_PRINT_LEN, PERC_PRINT_NAME);
    }
    // Print custom calculations on screen and in info file
    if (nCustomValues > 0) {
        print_customcalc(nPrintscr, printscrindex, customValues, customNames, MAX_PRINT_LEN, PERC_PRINT_NAME);
        fprint_customcalc(info, nPrintscr, printscrindex, customValues, customNames, MAX_PRINT_LEN, PERC_PRINT_NAME);
    }    
}