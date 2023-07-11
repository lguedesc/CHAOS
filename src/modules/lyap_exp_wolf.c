#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../libs/odesystems.h"
#include "../libs/nldyn.h"
#include "../libs/iofiles.h"
#include "../libs/interface.h"
#include "../libs/defines.h"
#include "../libs/basic.h"
#include "lyap_exp_wolf.h"

static void read_params(int dim, int npar, int *np, int *ndiv, int* trans, double *t, double **par, double **x);
static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double h, double t, double *x, double *par, char* funcname, size_t maxlength, double percname, char* mode);

void lyapunov_exp_wolf(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *)) {
    
    // Declare Program Parameters
    int nP;                         // Number of forcing periods analyzed
    int nDiv;                       // Number of divisions in each forcing period
    int trans;                      // Value of nP in which greater values are considered transient response
    // Assign values for program parameters, system parameters and initial conditions
    double t;
    double *x = NULL;
    double *par = NULL;
    read_params(DIM, nPar, &nP, &nDiv, &trans, &t, &par, &x);
    // Define Timestep
    double h = (2 * PI) / (nDiv * par[0]); // par[0] = OMEGA
    // Create output files to store results
    const char *directory = "data/LyapunovExp/out/";
    const char *module = "lyap";
    FILE *output_lyap = name_and_create_output_files(outputname, directory, module, ".csv");
    FILE *output_info = name_and_create_output_files(outputname, directory, "info", ".txt");    
    // Print information in screen and info output file
    print_info(output_info, DIM, nPar, nP, nDiv, trans, h, t, x, par, funcname, MAX_PRINT_LEN, PERC_PRINT_NAME, "screen");
    print_info(output_info, DIM, nPar, nP, nDiv, trans, h, t, x, par, funcname, MAX_PRINT_LEN, PERC_PRINT_NAME, "file");
    // Call solution
    lyap_wolf_solution(output_lyap, DIM, nP, nDiv, trans, t, &x, h, par, edosys, write_results_lyap);   
    // Close output file
    close_files(2, output_lyap, output_info);
    // Free allocated memory
    free_mem(x, par, NULL);
}

static void read_params(int dim, int npar, int *np, int *ndiv, int* trans, double *t, double **par, double **x) {
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
    // Close input file
    fclose(input);
    // Free memory
    free(input_filename);
    /* The user is responsible to free (x) and (par) after the function call */
}

static void print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double h, double t, double *x, double *par, char* funcname, size_t maxlength, double percname, char* mode) {
    //Get time and date
    time_t tm;
    time(&tm);

    if (strcmp(mode, "screen") == 0) {   
        write_prog_parameters_lyapunov(dim, npar, np, ndiv, trans, h, maxlength, percname);
        write_initial_conditions(dim, x, t, maxlength, percname);
        write_sys_parameters(npar, par, maxlength, percname);
        partition(2, maxlength); 
    } 
    else if (strcmp(mode, "file") == 0) {
        fwrite_prog_parameters_lyapunov(info, funcname, dim, npar, np, ndiv, trans, h, maxlength, percname);
        fwrite_initial_conditions(info, dim, x, t, maxlength, percname);
        fwrite_sys_parameters(info, npar, par, maxlength, percname);
        fpartition(info, 2, maxlength);
    }
    else {
        printf("Information could not be printed using mode (%s)...\n", mode);
        return;
    }

}
