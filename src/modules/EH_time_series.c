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
#include "EH_time_series.h"

void EH_timeseries(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *)) {
    
    // Declare Program Parameters
    const double pi = 4 * atan(1);  // Pi number definition
    int DIM;                        // Dimension of the system
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int nPar;                       // Number of parameters of the system
    int trans;                      // Number of forcing periods considered transient
    int nRMS;                       // Number of state variables that will be submitted to RMS calculation
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;                       // Time
    double *x = NULL;               // State Variables
    double *par = NULL;             // Parameters
    int *rmsindex = NULL;           // Indexes of state variables that will be submitted to RMS calculation
    double *xRMS = NULL;            // State Variables RMS values at permanent regime
    double *overallxRMS = NULL;     // State Variables RMS values at transient + permanent regime
    
    EH_tseries_read_params_and_IC(input_filename, &DIM, &nPar, &nP, &nDiv, &trans, &nRMS, &t, &par, &x, &rmsindex);
    
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
    EH_tseries_print_info(output_info, DIM, nPar, nP, nDiv, trans, nRMS, h, t, x, par, rmsindex, funcname,"screen");
    EH_tseries_print_info(output_info, DIM, nPar, nP, nDiv, trans, nRMS, h, t, x, par, rmsindex, funcname, "file");
    
    /*
    // Time variables
    double time_spent = 0.0;
    clock_t time_i = clock();
    */
    // Call solution
    EH_rk4_solution(output_rk4, DIM, nP, nDiv, trans, t, x, h, par, nRMS, rmsindex, &xRMS, &overallxRMS, edosys, write_results);

    /*
    clock_t time_f = clock();
    time_spent += (double)(time_f - time_i) / CLOCKS_PER_SEC; 
    printf("The elapsed time is %f seconds\n", time_spent);
    */

   // Print RMS calculations in screen and in info file
    EH_print_RMS(output_info, nRMS, rmsindex, xRMS, overallxRMS);
    
    // Close output file
    fclose(output_rk4);
    fclose(output_info);
    
    // Free allocated memory
    free(dir);
    free(input_filename);
    free(x); free(par);
    free(xRMS); free(overallxRMS);
}

void EH_tseries_read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, int *nrms, double *t, double **par, double **x, int **rmsindex) {
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
    // Close input file
    fclose(input);
    /* The user is responsible to free (x) and (par) after the function call */
}

void EH_tseries_print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, int nrms, double h, double t, double *x, double *par, int *rmsindex, char* funcname, char* mode) {
    //Get time and date
    time_t tm;
    time(&tm);

    if (strcmp(mode, "screen") == 0) {   
        printf("\n  Program Parameters\n");
        printf("  -------------------------------------------------\n");
        printf("%-30s%s%-20d\n", "  Dimension:", " ", dim);
        printf("%-30s%s%-20d\n", "  Number of Parameters:", " ", npar);
        printf("%-30s%s%-20d\n", "  Forcing Periods:", " ", np);
        printf("%-30s%s%-20d\n", "  Timesteps per Period:", " ", ndiv);
        printf("%-30s%s%-20d\n", "  Transient Considered:", " ", trans);
        printf("%-30s%s%-20g\n", "  Timestep value:", " ", h);
        printf("  -------------------------------------------------\n");
        printf("  Initial Conditions\n");
        printf("  -------------------------------------------------\n");
        printf("%-30s%s%-20g\n", "  Initial Time (t):", " ",  t);
        for (int i = 0; i < dim; i++) {
            printf("%s%d%-25s%s%-20g\n", "  x[", i, "]:", " ", x[i]);
        }
        printf("  -------------------------------------------------\n");
        printf("  System Parameters\n");
        printf("  -------------------------------------------------\n");
        for (int i = 0; i < npar; i++) {
            printf("%s%d%-23s%s%-20g\n", "  par[", i, "]:", " ", par[i]);
        }
        printf("  -------------------------------------------------\n");
        printf("  RMS Calculation Parameters\n");
        printf("  -------------------------------------------------\n");
        printf("%-30s%s%-20d\n", "  Number of RMS calculations:", " ", nrms);
        for (int i = 0; i < nrms; i++) {
            printf("%s%d%-18s%s%-20d\n", "  rmsindex[", i, "]:", " ", rmsindex[i]);
        }
    } 
    else if (strcmp(mode, "file") == 0) {
        fprintf(info, "  Date/Time:  %s", ctime(&tm)); 
        fprintf(info, "\n  =================================================\n");
        fprintf(info, "  Time Series: %s\n", funcname);
        fprintf(info, "  =================================================\n\n");
        fprintf(info, "  Program Parameters\n");
        fprintf(info, "  -------------------------------------------------\n");
        fprintf(info, "%-30s%s%-20d\n", "  Dimension:", " ", dim);
        fprintf(info, "%-30s%s%-20d\n", "  Number of Parameters:", " ", npar);
        fprintf(info, "%-30s%s%-20d\n", "  Forcing Periods:", " ", np);
        fprintf(info, "%-30s%s%-20d\n", "  Timesteps per Period:", " ", ndiv);
        fprintf(info, "%-30s%s%-20d\n", "  Transient Considered:", " ", trans);
        fprintf(info, "%-30s%s%-20g\n", "  Timestep value:", " ", h);
        fprintf(info, "  -------------------------------------------------\n");
        fprintf(info, "  Initial Conditions\n");
        fprintf(info, "  -------------------------------------------------\n");
        fprintf(info, "%-30s%s%-20g\n", "  Initial Time (t):", " ",  t);
        for (int i = 0; i < dim; i++) {
            fprintf(info, "%s%d%-25s%s%-20g\n", "  x[", i, "]:", " ", x[i]);
        }
        fprintf(info, "  -------------------------------------------------\n");
        fprintf(info, "  System Parameters\n");
        fprintf(info, "  -------------------------------------------------\n");
        for (int i = 0; i < npar; i++) {
            fprintf(info, "%s%d%-23s%s%-20g\n", "  par[", i, "]:", " ", par[i]);
        }
        fprintf(info, "  -------------------------------------------------\n");
        fprintf(info, "  RMS Calculation Parameters\n");
        fprintf(info, "  -------------------------------------------------\n");
        fprintf(info, "%-30s%s%-20d\n", "  Number of RMS calculations:", " ", nrms);
        for (int i = 0; i < nrms; i++) {
            fprintf(info, "%s%d%-18s%s%-20d\n", "  rmsindex[", i, "]:", " ", rmsindex[i]);
        }

    }
    else {
        printf("Information could not be printed using mode (%s)...\n", mode);
        return;
    }
}
