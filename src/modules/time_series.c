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
#include "time_series.h"

static void read_params(int dim, int npar, int *np, int *ndiv, double *t, double **par, double **x);
static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, double h, double t, double *x, double *par, char* funcname, char* mode);

void timeseries(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *)) {
    
    // Declare Program Parameters
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    // Assign values for program parameters, system parameters and initial conditions
    double t;
    double *x = NULL;
    double *par = NULL;
    read_params(DIM, nPar, &nP, &nDiv, &t, &par, &x);
    // Define Timestep
    double h = (2 * PI) / (nDiv * par[0]); // par[0] = OMEGA
    // Create output files to store results
    const char *directory = "data/TimeSeries/out/";                                            // Directory of output file
    const char *module = "timeseries";
    FILE *output_timeseries = name_and_create_output_files(outputname, directory, module, ".csv");    // Create timeseries output file 
    FILE *output_info = name_and_create_output_files(outputname, directory, "info", ".txt");   // Create info output file
    // Print information in screen and info output file
    print_info(output_info, DIM, nPar, nP, nDiv, h, t, x, par, funcname, "screen");
    print_info(output_info, DIM, nPar, nP, nDiv, h, t, x, par, funcname, "file");
    // Call solution
    rk4_solution(output_timeseries, DIM, nP, nDiv, t, x, h, par, edosys, write_results);
    // Close output file
    close_files(2, output_timeseries, output_info);        
    // Free allocated memory
    free_mem(x, par, NULL);
}

static void read_params(int dim, int npar, int *np, int *ndiv, double *t, double **par, double **x) {
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

static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, double h, double t, double *x, double *par, char* funcname, char* mode) {
    //Get time and date
    time_t tm;
    time(&tm);
    size_t maxlen = 100;

    if (strcmp(mode, "screen") == 0) {   
        printf("\n  Program Parameters\n");
        printf("  -------------------------------------------------\n");
        printf("%-30s%s%-20d\n", "  Dimension:", " ", dim);
        printf("%-30s%s%-20d\n", "  Number of Parameters:", " ", npar);
        printf("%-30s%s%-20d\n", "  Forcing Periods:", " ", np);
        printf("%-30s%s%-20d\n", "  Timesteps per Period:", " ", ndiv);
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
    }
    else {
        printf("Information could not be printed using mode (%s)...\n", mode);
        return;
    }

}
