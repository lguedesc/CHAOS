#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "bifurcation.h"

static void read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *parindex, double **x, int *bifmode);
static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int parindex, int bifmode, char* funcname, char* mode);

void bifurcation(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *)) {
    // Declare Program Parameters
    const double pi = 4 * atan(1);  // Pi number definition
    int DIM;                        // Dimension of the system
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int nPar;                       // Number of parameters of the system
    int trans;                      // Value of nP in which greater values are considered transient response
    int bMode;                      // Mode of the Bifurcation Diagram: 0 to follow attractor, 1 to reset initial conditions in each step
    int parIndex;                   // Index of control parameter in par array
    // Assign values for program parameters, system parameters and initial conditions
    char *input_filename = get_input_filename();
    double t;
    double *x = NULL;
    double *par = NULL;
    double *parRange = NULL;
    read_params_and_IC(input_filename, &DIM, &nPar, &nP, &nDiv, &trans, &t, &par, &parRange, &parIndex, &x, &bMode);
    
    // Create output files to store results
    char output_bifurc_name[200];
    char output_info_name[200];
    char output_poinc_name[200];
    const char *rawdir = "data/Bifurcation/out/";                                                            // Directory of output file
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
    print_info(output_info, DIM, nPar, nP, nDiv, trans, t, x, par, parRange, parIndex, bMode, funcname, "screen");
    print_info(output_info, DIM, nPar, nP, nDiv, trans, t, x, par, parRange, parIndex, bMode, funcname, "file");
    // To store the runtime of the program
    //double time_spent = 0.0;
    //clock_t time_i = clock();
    // Call solution
    bifurc_solution(output_bifurc, output_bifurc_poinc, DIM, nP, nDiv, trans, t, x, parIndex, parRange, par, edosys, write_bifurc_results, bMode);
    // Close output file
    fclose(output_bifurc);
    fclose(output_info);

    // Free allocated memory
    free(dir);
    free(input_filename);
    free(x); free(par); free(parRange);

    // Calculate time of execution
    //clock_t time_f = clock();
    //time_spent += (double)(time_f - time_i) / CLOCKS_PER_SEC; 
    //printf("The elapsed time is %f seconds", time_spent);

}

static void read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *parindex, double **x, int *bifmode) {
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
    // Read and assign program parameters
    fscanf(input, "%i %i %i", np, ndiv, trans); 
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
    // Determine the mode of the bifurcation diagram: 0 to follow attractor, 1 to reset initial conditions in each step
    fscanf(input, "%i", &(*bifmode));
    // Close input file
    fclose(input);
    /* The user is responsible to free (x), (par) or (parrange) after the function call */
}

static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int parindex, int bifmode, char* funcname, char* mode) {
    //Get time and date
    time_t tm;
    time(&tm);

    if (strcmp(mode, "screen") == 0) {   
        if (bifmode == 0){
            printf("\n  Bifurcation Mode: Following Attractor\n");
        }
        else if (bifmode == 1) {
            printf("\n  Bifurcation Mode: Reseting ICs\n");
        }
        else {
            printf("\n  Invalid Bifurcation Mode of %d...\nCheck Results!\n", bifmode);
        }
        printf("  -------------------------------------------------\n");
        printf("  Program Parameters\n");
        printf("  -------------------------------------------------\n");
        printf("%-30s%s%-20d\n", "  Dimension:", " ", dim);
        printf("%-30s%s%-20d\n", "  Number of Parameters:", " ", npar);
        printf("%-30s%s%-20d\n", "  Forcing Periods:", " ", np);
        printf("%-30s%s%-20d\n", "  Timesteps per Period:", " ", ndiv);
        printf("%-30s%s%-20d\n", "  Transient Considered:", " ", trans);
        printf("%-30s%s%-20s\n", "  Timestep value:", " ", "(2*pi)/(nDiv*par[0])");
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
        printf("%-30s%s%-20d\n", "  Parameter Index:", " ", parindex);
        printf("%-30s%s%-20g\n", "  Intial Parameter:", " ", parrange[0]);
        printf("%-30s%s%-20g\n", "  Final Parameter:", " ", parrange[1]);
        printf("%-30s%s%-20g\n", "  Increment Parameter:", " ", (parrange[1] - parrange[0]) / (parrange[2] - 1));
        printf("%-30s%s%-20g\n", "  Number of Steps:", " ", parrange[2]);
        printf("  -------------------------------------------------\n");        
    } 
    else if (strcmp(mode, "file") == 0) {
        fprintf(info, "  Date/Time:  %s", ctime(&tm)); 
        fprintf(info, "\n  ===================================================\n");
        fprintf(info, "  Bifurcation Diagram: %s\n", funcname);
        if (bifmode == 0){
            fprintf(info, "  Bifurcation Mode: Following Attractor\n");
        }
        else if (bifmode == 1) {
            fprintf(info, "  Bifurcation Mode: Reseting ICs\n");
        }
        else {
            fprintf(info, "  Invalid Bifurcation Mode of %i...\nCheck Results!\n", bifmode);
        }
        fprintf(info, "  ===================================================\n\n");
        fprintf(info, "\n  Program Parameters\n");
        fprintf(info, "  -------------------------------------------------\n");
        fprintf(info, "%-30s%s%-20d\n", "  Dimension:", " ", dim);
        fprintf(info, "%-30s%s%-20d\n", "  Number of Parameters:", " ", npar);
        fprintf(info, "%-30s%s%-20d\n", "  Forcing Periods:", " ", np);
        fprintf(info, "%-30s%s%-20d\n", "  Timesteps per Period:", " ", ndiv);
        fprintf(info, "%-30s%s%-20d\n", "  Transient Considered:", " ", trans);
        fprintf(info, "%-30s%s%-20s\n", "  Timestep value:", " ", "(2*pi)/(nDiv*par[0])");
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
        fprintf(info, "%-30s%s%-20d\n", "  Parameter Index:", " ", parindex);
        fprintf(info, "%-30s%s%-20g\n", "  Intial Parameter:", " ", parrange[0]);
        fprintf(info, "%-30s%s%-20g\n", "  Final Parameter:", " ", parrange[1]);
        fprintf(info, "%-30s%s%-20g\n", "  Increment Parameter:", " ", (parrange[1] - parrange[0]) / (parrange[2] - 1));
        fprintf(info, "%-30s%s%-20g\n", "  Number of Steps:", " ", parrange[2]);
        fprintf(info, "  -------------------------------------------------\n");
    }
    else {
        printf("Information could not be printed using mode (%s)...\n", mode);
        return;
    }

}
