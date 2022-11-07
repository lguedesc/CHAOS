#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../libs/interface.h"
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "../libs/energyharvest.h"
#include "EH_time_series.h"

static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, int nrms, double h, double t, double *x, double *par, int *rmsindex, char* funcname, 
                       int ncustomvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex, size_t maxlength, double percname, char* mode);
static void read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, int *nrms, double *t, double **par, double **x, int **rmsindex, int *ncustomvalues, int *nprintf, int *nprintscr, int **printfindex, int **printscrindex);

void EH_timeseries(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, int, int, char **, size_t, double *, int)) {
    
    // Parameters related to printing information
    size_t maxLen = 71;             // Max length of the info printed on the screen and on info file
    double percName = 0.6;          // Percentage of space occuped by the name of the quantity printed
    // Declare Program Parameters related to the system of equations
    const double pi = 4 * atan(1);  // Pi number definition
    int DIM;                        // Dimension of the system
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int nPar;                       // Number of parameters of the system
    int trans;                      // Number of forcing periods considered transient
    int nRMS;                       // Number of state variables that will be submitted to RMS calculation
    int nCustomValues = 0;          // If there is a custom function to be called, this is the number of calculations the function is going to perform
    int nPrintf = 0;                // Number of custom values to be printed in the output file
    int nPrintscr = 0;              // Number of custom values to be printed on the screen
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;                           // Time
    double *x = NULL;                   // State Variables
    double *par = NULL;                 // Parameters
    int *rmsindex = NULL;               // Indexes of state variables that will be submitted to RMS calculation
    double *xRMS = NULL;                // State Variables RMS values at permanent regime
    double *overallxRMS = NULL;         // State Variables RMS values at transient + permanent regime
    double *customValues = NULL;        // Variable to store all custom values that custom functions calculate
    char **customnames = NULL;          // Names of the custom values
    int *printfindex = NULL;            // Indexes of custom values that will be printed in the output file
    int *printscrindex = NULL;          // Indexes of custom values that will be printed on the screen

    read_params_and_IC(input_filename, &DIM, &nPar, &nP, &nDiv, &trans, &nRMS, &t, &par, &x, &rmsindex, &nCustomValues, &nPrintf, &nPrintscr, &printfindex, &printscrindex);
    // Define Timestep
    double h = (2 * pi) / (nDiv * par[0]); // par[0] = OMEGA
    
    // Create output files to store results
    char output_rk4_name[200];
    char output_info_name[200];
    const char *rawdir = "TimeSeries/out/";                                                              // Directory of output file
    char *dir = convert_dir(rawdir);
    const char *ext = ".csv";                                                                           // Extension of output file    
    const char *ext_info = ".txt";                                                                      // Extension of info file
    snprintf(output_rk4_name, sizeof(output_rk4_name), "%s%s_rk4", dir, outputname);                      // Assign name for output rk4 without extension
    snprintf(output_info_name, sizeof(output_info_name), "%s%s_info", dir, outputname);                   // Assign name for output info without extension
    FILE *output_rk4 = create_output_file(output_rk4_name, ext, dir);    // Create rk4 output file 
    FILE *output_info = create_output_file(output_info_name, ext_info, dir);  // Create info output file
    
    // Print information in screen and info output file
    print_info(output_info, DIM, nPar, nP, nDiv, trans, nRMS, h, t, x, par, rmsindex, funcname, nCustomValues, nPrintf, printfindex, nPrintscr, printscrindex, maxLen, percName, "screen");
    print_info(output_info, DIM, nPar, nP, nDiv, trans, nRMS, h, t, x, par, rmsindex, funcname, nCustomValues, nPrintf, printfindex, nPrintscr, printscrindex, maxLen, percName,"file");
    /*
    // Time variables
    double time_spent = 0.0;
    clock_t time_i = clock();
    */
    // Call solution
    
    EH_rk4_solution(output_rk4, DIM, nP, nDiv, trans, t, x, h, par, nRMS, rmsindex, &xRMS, &overallxRMS, edosys, EH_write_timeseries_results, 
                    nCustomValues, &customnames, &customValues, nPrintf, printfindex, nPrintscr, printscrindex, customfunc);
    /*
    clock_t time_f = clock();
    time_spent += (double)(time_f - time_i) / CLOCKS_PER_SEC; 
    printf("The elapsed time is %f seconds\n", time_spent);
    */
    
    // Print RMS calculations in screen and in info file
    print_RMS(nRMS, rmsindex, xRMS, overallxRMS, maxLen, percName);
    fprint_RMS(output_info, nRMS, rmsindex, xRMS, overallxRMS, maxLen, percName);
    // Print custom calculations on screen and in info file
    for (int i = 0; i < nCustomValues; i++) {
        printf("customnames[%d] = %s\n", i, customnames[i]);
    }
    if (nCustomValues > 0) {
        print_customcalc(nPrintscr, printscrindex, customValues, customnames, maxLen, percName);
        fprint_customcalc(output_info, nPrintscr, printscrindex, customValues, customnames, maxLen, percName);
    }    
    // Close output file
    fclose(output_rk4);
    fclose(output_info);
    
    // Free allocated memory
    free(dir);
    free(input_filename);
    free(x); free(par);
    free(xRMS); free(overallxRMS); free(rmsindex);
    if (nCustomValues > 0) {
        free(customValues); free(printfindex); free(printscrindex);
    }
}

static void read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, int *nrms, double *t, double **par, double **x, int **rmsindex, int *ncustomvalues, int *nprintf, int *nprintscr, int **printfindex, int **printscrindex) {
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
