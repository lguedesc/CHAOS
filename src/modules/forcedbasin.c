#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/edosystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "forcedbasin.h"

void forcedbasin(char *funcname, void (*edosys)(int, double *, double, double *, double *)) {
    // Declare Program Parameters
    const double pi = 4 * atan(1);  // Pi number definition
    int DIM;                        // Dimension of the system
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int nPar;                       // Number of parameters of the system
    int trans;                      // Value of nP in which greater values are considered transient response
    int maxPer;                     // Maximum periodicity to be classified
    int indexX;                     // Index of control parameter in X
    int indexY;                     // Index of control parameter in Y
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;
    double *x = NULL;
    double *par = NULL;
    double *icRange = NULL;
    forcedbasin_read_params_and_IC(input_filename, &DIM, &nPar, &maxPer, &nP, &nDiv, &trans, &t, &par, &icRange, &indexX, &indexY, &x);
    
    // Create output files to store results
    char output_dyndiag_name[200];
    char output_info_name[200];
    const char *rawdir = "ForcBasin/out/";                                                            // Directory of output file
    char *dir = convert_dir(rawdir);
    const char *ext = ".csv";                                                                           // Extension of output file
    const char *ext_info = ".txt";                                                                      // Extension of info file
    snprintf(output_dyndiag_name, sizeof(output_dyndiag_name), "%s%s_forcedbasin", dir, funcname);            // Assign name for output file without extension
    snprintf(output_info_name, sizeof(output_info_name), "%s%s_info", dir, funcname);                   // Assign name for output info without extension
    FILE *output_dyndiag = create_output_file(output_dyndiag_name, ext, dir);                             // Create dynamical diagram output file 
    FILE *output_info = create_output_file(output_info_name, ext_info, dir);                            // Create info output file
    
    // Print information in screen and info output file
    forcedbasin_print_info(output_info, DIM, nPar, maxPer, nP, nDiv, t, x, par, icRange, indexX, indexY, funcname, "screen");
    forcedbasin_print_info(output_info, DIM, nPar, maxPer, nP, nDiv, t, x, par, icRange, indexX, indexY, funcname, "file");
    // To store the runtime of the program
    //double time_spent = 0.0;
    //clock_t time_i = clock();
    // Call solution
    forced_basin_of_attraction_2D(output_dyndiag, DIM, nP, nDiv, trans, maxPer, t, &x, indexX, indexY, icRange, par, nPar, edosys, p_write_dyndiag_results);
    
    // Close output file
    fclose(output_dyndiag);
    fclose(output_info);

    // Free allocated memory
    free(dir);
    free(input_filename);
    free(x); free(par); free(icRange);

    // Calculate time of execution
    //clock_t time_f = clock();
    //time_spent += (double)(time_f - time_i) / CLOCKS_PER_SEC; 
    //printf("The elapsed time is %f seconds", time_spent);

}

void forcedbasin_read_params_and_IC(char *name, int *dim, int *npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **icrange, int *indexX, int *indexY, double **x) {
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

void forcedbasin_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, double t, double *x, double *par, double *icrange, int indexX, int indexY, char* edosys, char* mode) {
    //Get time and date
    time_t tm;
    time(&tm);

    if (strcmp(mode, "screen") == 0) {   
        printf("\n==================================================\n");
        printf("%-30s%-20s\n", "Basin of Attraction (Forced): ", edosys);
        printf("%-30s%-4.0lf%-3s%-4.0lf\n", "Resolution: ", icrange[2], " x ", icrange[5]);
        printf("==================================================\n\n");
        printf("Program Parameters\n");
        printf("--------------------------------------------------\n");
        printf("%-30s%-12d\n", "Dimension: ", dim);
        printf("%-30s%-12d\n", "Number of Parameters: ", npar);
        printf("%-30s%-12d\n", "Max Periodicity Class: ", maxper);
        printf("%-30s%-12d\n", "Forcing Periods: ", np);
        printf("%-30s%-12d\n", "Timesteps per Period: ", ndiv);
        printf("%-30s%-20s\n", "Timestep value: ", "(2*pi)/(nDiv*par[0])");
        printf("--------------------------------------------------\n");
        printf("Initial Conditions\n");
        printf("--------------------------------------------------\n");
        printf("%-30s%-12.10lf\n", "Initial Time (t): ", t);
        for (int i = 0; i < dim; i++) {
            printf("x[%d]:%-25s%.10lf\n", i, "", x[i]);
        }
        printf("--------------------------------------------------\n");
        printf("System Parameters\n");
        printf("--------------------------------------------------\n");
        for (int i = 0; i < npar; i++) {
            printf("par[%d]:%-23s%-12g\n", i, "", par[i]);
        }
        printf("--------------------------------------------------\n");
        printf("%-30s%-12d\n", "Parameter Index (x): ", indexX);
        printf("%-30s%-12g\n", "Initial Parameter (x): ", icrange[0]);
        printf("%-30s%-12g\n", "Final Parameter (x): ", icrange[1]);
        printf("%-30s%-12g\n", "Increment Parameter (x): ", (icrange[1] - icrange[0])/(icrange[2] - 1));
        printf("--------------------------------------------------\n");
        printf("%-30s%-12d\n", "Parameter Index (y): ", indexY);
        printf("%-30s%-12g\n", "Initial Parameter (y): ", icrange[3]);
        printf("%-30s%-12g\n", "Final Parameter (y): ", icrange[4]);
        printf("%-30s%-12g\n", "Increment Parameter (y): ", (icrange[4] - icrange[3])/(icrange[5] - 1));
        printf("--------------------------------------------------\n");
        
    } 
    else if (strcmp(mode, "file") == 0) {
        fprintf(info, "Date/Time:  %s", ctime(&tm)); 
        fprintf(info, "\n==================================================\n");
        fprintf(info, "%-30s%-20s\n", "Basin of Attraction (Forced): ", edosys);
        fprintf(info, "%-30s%-4.0lf%-3s%-4.0lf\n", "Resolution: ", icrange[2], " x ", icrange[5]);
        fprintf(info, "==================================================\n\n");
        fprintf(info, "Program Parameters\n");
        fprintf(info, "--------------------------------------------------\n");
        fprintf(info, "%-30s%-20d\n", "Dimension: ", dim);
        fprintf(info, "%-30s%-20d\n", "Number of Parameters: ", npar);
        fprintf(info, "%-30s%-20d\n", "Max Periodicity Class: ", maxper);
        fprintf(info, "%-30s%-20d\n", "Forcing Periods: ", np);
        fprintf(info, "%-30s%-20d\n", "Timesteps per Period: ", ndiv);
        fprintf(info, "%-30s%-20s\n", "Timestep value: ", "(2*pi)/(nDiv*par[0])");
        fprintf(info, "--------------------------------------------------\n");
        fprintf(info, "Initial Conditions\n");
        fprintf(info, "--------------------------------------------------\n");
        fprintf(info, "%-30s%-20.10lf\n", "Initial Time (t): ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(info, "x[%d]:%-25s%.10lf\n", i, "", x[i]);
        }
        fprintf(info, "--------------------------------------------------\n");
        fprintf(info, "System Parameters\n");
        fprintf(info, "--------------------------------------------------\n");
        for (int i = 0; i < npar; i++) {
            fprintf(info, "par[%d]:%-23s%-12g\n", i, "", par[i]);
        }
        fprintf(info, "--------------------------------------------------\n");
        fprintf(info, "%-30s%-12d\n", "Parameter Index (x): ", indexX);
        fprintf(info, "%-30s%-12g\n", "Initial Parameter (x): ", icrange[0]);
        fprintf(info, "%-30s%-12g\n", "Final Parameter (x): ", icrange[1]);
        fprintf(info, "%-30s%-12g\n", "Increment Parameter (x): ", (icrange[1] - icrange[0])/(icrange[2] - 1));
        fprintf(info, "--------------------------------------------------\n");
        fprintf(info, "%-30s%-12d\n", "Parameter Index (y): ", indexY);
        fprintf(info, "%-30s%-12g\n", "Initial Parameter (y): ", icrange[3]);
        fprintf(info, "%-30s%-12g\n", "Final Parameter (y): ", icrange[4]);
        fprintf(info, "%-30s%-12g\n", "Increment Parameter (y): ", (icrange[4] - icrange[3])/(icrange[5] - 1));
        fprintf(info, "--------------------------------------------------\n");
    }
    else {
        printf("Information could not be printed in file using mode (%s)...\n", mode);
        return;
    }

}
