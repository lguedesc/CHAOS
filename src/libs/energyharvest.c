#include <stdio.h>
#include <stdlib.h>
#include "edosystems.h"
#include "iofiles.h"
#include "nldyn.h"

double dimensionless_inst_out_power(double chi, double varphi, double kappa, double v) {
    return (((chi * varphi)/kappa)*v*v);
}

double dimensionless_RMS(double *cum, double measure, int N, int mode) {
    if (mode == 0) {
        // accumulate the value of the square of the measure 
        (*cum) = (*cum) + (measure * measure);
        return (*cum);
    }
    else if (mode == 1) {
        double RMS;
        // Computes the RMS value (cummulative / N times accumulated)
        RMS = sqrt( (*cum) / N );
        return RMS;
    }
    else {
        printf("Failed to compute dimensionless RMS using mode (%i)\n", mode);
        return;
    }
}

void write_results_EH_ftimeseries(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, double Pout, int mode) {
    // Check the mode of the function
    if (mode == 1) {
        // Header
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%i] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "sLE[%i] ", i);
        }
        fprintf(output_file, "Pout\n");
        // Initial Conditions
        fprintf(output_file, "%.10lf ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", s_lambda[i]);
        }
        fprintf(output_file, "%.10lf\n", Pout);
    }
    else if (mode == 2) {
        fprintf(output_file, "%.10lf ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", s_lambda[i]);
        }
        fprintf(output_file, "%.10lf\n", Pout);
    } 
    else if (mode == 3) {
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 4) {
        fprintf(output_file, "%.10lf ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        fprintf(output_file, "\n");
    }
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void EH_full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, double Pout, int mode)) {
    // Allocate memory for x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Allocate memory for vectors necessary for lyapunov exponents calculation
    double *cum = malloc(dim * sizeof *cum);            // Cumulative Vector
    double *lambda = malloc(dim *sizeof *lambda);       // Lyapunov Exponents vector
    double *s_cum = malloc(dim * sizeof *s_cum);        // Short Cumulative Vector
    double *s_lambda = malloc(dim * sizeof *s_lambda);  // Short Lyapunov Exponents Vector
    double *znorm = malloc(dim * sizeof *znorm);        // Norm of Vectors
    double *gsc = malloc((dim - 1) * sizeof *gsc);      // Inner Products Vector
    // Assign initial values to Lyapunov exponents vector
    for (int i = 0; i < dim; i++) {
        lambda[i] = 0.0;
        s_lambda[i] = 0.0;
    }
    // Numerical control parameters
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations
    int N = dim * ndiv;                                 // Number of steps in integration
    double tf = h*np*ndiv;                              // Final time
    double s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time 
    // Declare vector and allocate memory to store poincare map values: poinc[number of permanent regime forcing periods][dimension original system]
    double **poinc = malloc((np - trans) * sizeof **poinc);
    for (int i = 0; i < np - trans; i++) {
        poinc[i] = malloc(dim * sizeof **poinc);
    }
    // Declare vector to store the chosen Lyapunov Exponents to determine the attractor
    double *LE = malloc(dim * sizeof *LE);
    // Prepare x vector to include perturbed values 
    realloc_vector(x, ndim);
    // Assign initial perturbation
    perturb_wolf(x, dim, ndim, &cum, &s_cum);
    printf("after perturb (outside function)\n");
    for(int i = 0; i < ndim; i++) {
        printf("x[%i] = %lf (memory: %p)\n", i, (*x)[i], &(*x)[i]);
    }
    // Declare variables to calculate RMS output Power
    double Pout = 0, Pout_RMS = 0;
    // Make the header of output files
    write_results(output_ftimeseries_file, dim, t, (*x), lambda, s_lambda, Pout, 1);
    write_results(output_poinc_file, dim, t, (*x), lambda, s_lambda, Pout, 3);
    // Call Runge-Kutta 4th order integrator N = nP * nDiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(ndim, *x, t, h, par, f, edosys);
            lyapunov_wolf(x, t, h, dim, ndim, s_T0, &cum, &s_cum, &lambda, &s_lambda, &znorm, &gsc);
            t = t + h;
            // Apply poincare map at permanent regime
            if (i >= trans) {
                // Choose any point in the trajectory for poincare section placement
                if (j == 1) {
                    // Stores poincare values in poinc[np - trans][dim] vector
                    for (int p = 0; p < dim; p++) {
                        poinc[i - trans][p] = (*x)[p];
                    }
                    write_results(output_poinc_file, dim, t, (*x), lambda, s_lambda, Pout, 4);
                }
            // Compute Instantaneous Power and Squared values for RMS
            Pout = dimensionless_inst_out_power(par[5], par[6], par[7], (*x)[2]);
            Pout_RMS = dimensionless_RMS(&Pout_RMS, Pout, N, 0);
            }
            write_results(output_ftimeseries_file, dim, t, (*x), lambda, s_lambda, Pout, 2);
        }
    }
    // Compute Final RMS output Power
    Pout_RMS = dimensionless_RMS(&Pout_RMS, Pout, N, 1);
    printf("Pout_RMS = %.10lf\n", Pout_RMS);
    // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
    store_LE(dim, lambda, s_lambda, LE);
    for (int i = 0; i< dim; i++) {
        printf("LE[%i] = %.10lf\n", i, LE[i]);
    }
    // Declare vector for temporary storage of periodicity values to check if all directions are equal
    int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
    // Declare variable to flag if all directions present same periodicity or not (0 = all the same, 1 = not the same)
    int diffAttrac = -1;
    // Verify the type of motion of the system
    (*attrac) = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
    printf("diffAttrac = %i\n", diffAttrac);
    // Free memory    
    free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
    free(znorm); free(gsc); free(LE); free(tmp_attrac);
    for (int i = 0; i < np - trans; i++) {
        free(poinc[i]);
    }
    free(poinc);
}