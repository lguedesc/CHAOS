#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/edosystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "ftime_series.h"

void ftime_series(char *funcname, void (*edosys)(int, double *, double, double *, double *)) {
    
    // Declare Program Parameters
    const double pi = 4 * atan(1);  // Pi number definition
    int DIM;                        // Dimension of the system
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int nPar;                       // Number of parameters of the system
    int trans;                      // Value of nP in which greater values are considered transient response
    int attractor;                  // Value of type of motion of the system
    int maxPer;                     // Maximum periodicity to be classified
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;
    double *x = NULL;
    double *par = NULL;
    fts_read_params_and_IC(input_filename, &DIM, &nPar, &maxPer, &nP, &nDiv, &trans, &t, &par, &x);
    // Define Timestep
    double h = (2 * pi) / (nDiv * par[0]); // par[0] = OMEGA
    // Create output files to store results
    char output_ftimeseries_name[200];
    char output_poinc_name[200];
    char output_info_name[200];
    const char *rawdir = "FTimeSeries/out/";                                                                // Directory of output file
    char *dir = convert_dir(rawdir);
    const char *ext = ".csv";                                                                               // Extension of output file
    const char *ext_info = ".txt";                                                                          // Extension of info file
    snprintf(output_ftimeseries_name, sizeof(output_ftimeseries_name), "%s%s_ftimeseries", dir, funcname);  // Assign name for output ftimeseries file without extension
    snprintf(output_poinc_name, sizeof(output_poinc_name), "%s%s_poinc", dir, funcname);                    // Assign name for output poincare file without extension
    snprintf(output_info_name, sizeof(output_info_name), "%s%s_info", dir, funcname);                       // Assign name for output info without extension
    FILE *output_ftimeseries = create_output_file(output_ftimeseries_name, ext, dir);                       // Create ftimeseries output file 
    FILE *output_poinc = create_output_file(output_poinc_name, ext, dir);                                   // Create poincare output file 
    FILE *output_info = create_output_file(output_info_name, ext_info, dir);                                     // Create info output file
    
    // Print information in screen and in info file
    fts_print_info(output_info, DIM, nPar, maxPer, nP, nDiv, h, t, x, par, funcname, "screen");
    fts_print_info(output_info, DIM, nPar, maxPer, nP, nDiv, h, t, x, par, funcname, "file");
    // Call solution
    full_timeseries_solution(output_ftimeseries, output_poinc, DIM, nP, nDiv, trans, &attractor, maxPer, t, &x, h, par, edosys, write_results_ftimeseries);
    // Print attractor in screen and in info file
    fts_print_attractor(output_info, attractor, maxPer, "screen");
    fts_print_attractor(output_info, attractor, maxPer, "file");
    
    // Close output file
    fclose(output_ftimeseries);
    fclose(output_poinc);
    fclose(output_info);

    // Free allocated memory
    free(dir);
    free(input_filename);
    free(x); free(par);
}

void fts_read_params_and_IC(char *name, int *dim, int *npar, int *maxper, int *np, int *ndiv, int* trans, double *t, double **par, double **x) {
   // Open input file
    FILE *input = fopen(name, "r");
    if (input == NULL) {
        // Return error if input does not exist 
        perror(name);
        exit(1);
    }
    // Read and assign system constants
    fscanf(input, "%i", dim);
    fscanf(input, "%i", npar);
    fscanf(input, "%i", maxper);
    // Read and assign program parameters
    fscanf(input, "%i %i %i", np, ndiv, trans); 
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
    // Close input file
    fclose(input);
    /* The user is responsible to free (x) and (par) after the function call */
}

void fts_print_attractor(FILE* info, int attrac, int maxper, char *mode) {
    if (strcmp(mode, "file") == 0) {
        if (attrac < maxper) { 
            fprintf(info, "%-24s%-7s%-5i\n", "Type of Motion: ", "Period-", attrac); 
        }
        else if (attrac == maxper) {
            fprintf(info, "%-24s%-7s%-1i%-4s\n", "Type of Motion: ", "Period-", attrac, "+");
        }
        else if (attrac == maxper + 1) {
            fprintf(info, "%-24s%-13s\n", "Type of Motion: ", "Chaotic");
        }
        else if (attrac == maxper + 2) {
            fprintf(info, "%-24s%-13s\n", "Type of Motion: ", "Hyperchaotic");
        }
        else if (attrac == maxper + 3) {
            fprintf(info, "%-24s%-13s\n", "Type of Motion: ", "Undefined Periodicity");
        }
        else {
            fprintf(info, "%-24s%-13s\n", "Type of Motion: ", "Undefined (Escape)");
        }
    }
    else if (strcmp(mode, "screen") == 0) {
        if (attrac < maxper) { 
            printf("%-24s%-7s%-5i\n", "Type of Motion: ", "Period-", attrac); 
        }
        else if (attrac == maxper) {
            printf("%-24s%-7s%-1i%-4s\n", "Type of Motion: ", "Period-", attrac, "+");
        }
        else if (attrac == maxper + 1) {
            printf("%-24s%-13s\n", "Type of Motion: ", "Chaotic");
        }
        else if (attrac == maxper + 2) {
            printf("%-24s%-13s\n", "Type of Motion: ", "Hyperchaotic");
        }
        else if (attrac == maxper + 3) {
            printf("%-24s%-13s\n", "Type of Motion: ", "Undefined Periodicity");
        }
        else {
            printf("%-24s%-13s\n", "Type of Motion: ", "Undefined (Escape)");
        }
    }
    
}

void fts_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, double h, double t, double *x, double *par, char* edosys, char* mode) {
    //Get time and date
    time_t tm;
    time(&tm);
    
    if (strcmp(mode, "screen") == 0) {   
        printf("\n====================================\n");
        printf("Full Time series: %s\n", edosys);
        printf("====================================\n\n");
        printf("Program Parameters\n");
        printf("------------------------------------\n");
        printf("%-24s%-13d\n", "Dimension: ", dim);
        printf("%-24s%-13d\n", "Number of Parameters: ", npar);
        printf("%-24s%-13d\n", "Max Periodicity Class: ", maxper);
        printf("%-24s%-13d\n", "Forcing Periods: ", np);
        printf("%-24s%-13d\n", "Timesteps per Period: ", ndiv);
        printf("%-24s%-13g\n", "Timestep value: ", h);
        printf("------------------------------------\n");
        printf("Initial Conditions\n");
        printf("------------------------------------\n");
        printf("%-24s%-13lf\n", "Initial Time (t): ", t);
        for (int i = 0; i < dim; i++) {
            printf("%-2s%-1d%-21s%-13g\n", "x[", i, "]: ", x[i]);
        }
        printf("------------------------------------\n");
        printf("System Parameters\n");
        printf("------------------------------------\n");
        for (int i = 0; i < npar; i++) {
            printf("%-4s%-1d%-19s%-13g\n", "par[", i, "]: ", par[i]);
        }
        printf("------------------------------------\n");
    } 
    else if (strcmp(mode, "file") == 0) {
        fprintf(info, "Date/Time:  %s", ctime(&tm)); 
        fprintf(info, "\n====================================\n");
        fprintf(info, "Full Time series: %s\n", edosys);
        fprintf(info, "====================================\n\n");
        fprintf(info, "Program Parameters\n");
        fprintf(info, "------------------------------------\n");
        fprintf(info, "%-24s%-13d\n", "Dimension: ", dim);
        fprintf(info, "%-24s%-13d\n", "Number of Parameters: ", npar);
        fprintf(info, "%-24s%-13d\n", "Max Periodicity Class: ", maxper);
        fprintf(info, "%-24s%-13d\n", "Forcing Periods: ", np);
        fprintf(info, "%-24s%-13d\n", "Timesteps per Period: ", ndiv);
        fprintf(info, "%-24s%-13.10lf\n", "Timestep value: ", h);
        fprintf(info, "------------------------------------\n");
        fprintf(info, "Initial Conditions\n");
        fprintf(info, "------------------------------------\n");
        fprintf(info, "%-24s%-13.10lf\n", "Initial Time (t): ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(info, "%-2s%-1d%-21s%-13.10lf\n", "x[", i, "]: ", x[i]);
        }
        fprintf(info, "------------------------------------\n");
        fprintf(info, "System Parameters\n");
        fprintf(info, "------------------------------------\n");
        for (int i = 0; i < npar; i++) {
            fprintf(info, "%-4s%-1d%-19s%-13g\n", "par[", i, "]: ", par[i]);
        }
        fprintf(info, "------------------------------------\n");
    }
    else {
        printf("Information could not be printed using mode (%s)...\n", mode);
        return;
    }

}
