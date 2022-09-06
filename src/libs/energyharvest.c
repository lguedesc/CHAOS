#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iofiles.h"
#include "nldyn.h"

double RMS(double *cum, double measure, int N, int mode) {
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
        printf("DEBUG: Failed to compute dimensionless RMS using mode (%d)\n", mode);
        return (*cum);
    }
}

void EH_print_RMS(FILE *info, int nRMS, int *rmsindex, double *xRMS, double *overallxRMS) {
    // Print RMS calculations on screen and in file
    for (int q = 0; q < nRMS; q++) {
        printf("%s%d%-22s%s%-20g\n", "  xRMS[", rmsindex[q], "]:", " ", xRMS[rmsindex[q]]);
        fprintf(info, "%s%d%-22s%s%-20g\n", "  xRMS[", rmsindex[q], "]:", " ", xRMS[rmsindex[q]]);
    }
    for (int q = 0; q < nRMS; q++) {
        printf("%s%d%-14s%s%-20g\n", "  Overall xRMS[", rmsindex[q], "]:", " ", overallxRMS[rmsindex[q]]);
        fprintf(info, "%s%d%-14s%s%-20g\n", "  Overall xRMS[", rmsindex[q], "]:", " ", overallxRMS[rmsindex[q]]);
    }
    printf("  -------------------------------------------------\n");
    fprintf(info, "  -------------------------------------------------\n");
}

// Solutions
void EH_rk4_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, int nrms, int *rmsindex, double **xrms, double **overallxrms, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, int mode)) {
    // Allocate x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Allocate RMS variables
    (*xrms) = malloc(dim * sizeof **xrms);
    (*overallxrms) = malloc(dim * sizeof **overallxrms);
    // Initialize xRMS[dim] and overallxRMS[dim] vector
    for (int i = 0; i < dim; i ++) {
        (*xrms)[i] = 0.0;
        (*overallxrms)[i] = 0.0;
    }
    // Mumber of integration steps
    int N = np*ndiv;
    // Make the header o output file
    write_results(output_file, dim, t, x, 1);
    // Call Runge-Kutta 4th order integrator n = np * ndiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(dim, x, t, h, par, f, edosys);
            t = t + h;
            // Permanent Regime Calculations
            if (i >= trans) {
                // Accumulate squared values to RMS computation in permanent regime
                for (int q = 0; q < nrms; q++) {
                    (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], x[rmsindex[q]], N, 0);
                }
            }
            // Accumulate squared values to RMS computation for all time domain
            for (int q = 0; q < nrms; q++) {
                (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], x[rmsindex[q]], N, 0);
            }
            write_results(output_file, dim, t, x, 2);
        }
    }
    // Compute RMS values of state variables
    for (int q = 0; q < nrms; q++) {
        (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], x[rmsindex[q]], N, 1);
        (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], x[rmsindex[q]], N, 1);
    }
    // Free Memory
    free(f); 
}

void EH_full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, int nrms, int *rmsindex, double **xrms, double **overallxrms, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int mode)) {
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
    // Allocate RMS variables
    (*xrms) = malloc(dim * sizeof **xrms);
    (*overallxrms) = malloc(dim * sizeof **overallxrms);
    // Initialize xRMS[dim] and overallxRMS[dim] vector
    for (int i = 0; i < dim; i ++) {
        (*xrms)[i] = 0.0;
        (*overallxrms)[i] = 0.0;
    }
    // Mumber of integration steps
    int N = np*ndiv;
    // Numerical control parameters
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations
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
    // Make the header of output files
    write_results(output_ftimeseries_file, dim, t, (*x), lambda, s_lambda, 1);
    write_results(output_poinc_file, dim, t, (*x), lambda, s_lambda, 3);
    // Call Runge-Kutta 4th order integrator N = nP * nDiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(ndim, *x, t, h, par, f, edosys);
            lyapunov_wolf(x, t, h, dim, ndim, s_T0, &cum, &s_cum, &lambda, &s_lambda, &znorm, &gsc);
            t = t + h;
            write_results(output_ftimeseries_file, dim, t, (*x), lambda, s_lambda, 2);
            // Apply poincare map at permanent regime
            if (i >= trans) {
                // Choose any point in the trajectory for poincare section placement
                if (j == 1) {
                    // Stores poincare values in poinc[np - trans][dim] vector
                    for (int p = 0; p < dim; p++) {
                        poinc[i - trans][p] = (*x)[p];
                    }
                    write_results(output_poinc_file, dim, t, (*x), lambda, s_lambda, 4);
                }
                // Accumulate squared values to RMS computation in permanent regime
                for (int q = 0; q < nrms; q++) {
                    (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 0);
                }
            }
            // Accumulate squared values to RMS computation for all time domain
            for (int q = 0; q < nrms; q++) {
                (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 0);
            }
        }
    }
    // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
    store_LE(dim, lambda, s_lambda, LE);
    /*for (int i = 0; i< dim; i++) {
        printf("  LE[%i] = %.10lf\n", i, LE[i]);
    }*/
    // Declare vector for temporary storage of periodicity values to check if all directions are equal
    int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
    // Declare variable to flag if all directions present same periodicity or not (0 = all the same, 1 = not the same)
    int diffAttrac = -1;
    // Verify the type of motion of the system
    (*attrac) = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
    //printf("  diffAttrac = %i\n", diffAttrac);
    // Compute RMS values of state variables
    for (int q = 0; q < nrms; q++) {
        (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
        (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
    }
    // Free memory    
    free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
    free(znorm); free(gsc); free(LE); free(tmp_attrac);
    for (int i = 0; i < np - trans; i++) {
        free(poinc[i]);
    }
    free(poinc);
}

void EH_bifurc_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, double t, double *x, int parindex, 
                        double *parrange, double *par, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                        void (*write_results)(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, int nrms, int *rmsindex, double *xrms, double *overallxrms, int mode), int bifmode) {
    // Allocate x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Store Initial Conditions
    double t0 = t;
    double *IC = malloc(dim * sizeof *IC);
    for (int i = 0; i < dim; i++) {
        IC[i] = x[i];
    }
    // Declare memory to store min and max values
    double *xmax = malloc(dim * sizeof *xmax);
    double *xmin = malloc(dim * sizeof *xmin);
    // Allocate RMS variables
    double *xrms = malloc(dim * sizeof *xrms);
    double *overallxrms = malloc(dim * sizeof *overallxrms);
    // Mumber of integration steps
    int N = np*ndiv;
    // Declare timestep and pi
    double h;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control parameters
    double varstep;        
    varstep = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    // Make the header of the output file
    write_results(output_poinc_file, dim, par[parindex], x, xmin, xmax, nrms, rmsindex, xrms, overallxrms, 0);
    write_results(output_file, dim, par[parindex], x, xmin, xmax, nrms, rmsindex, xrms, overallxrms, 2);
    // Starts sweep the control parameter
    for (int k = 0; k < (int)parrange[2]; k++) {
        par[parindex] = parrange[0] + k*varstep; // Increment value
        // Reset Variables
        t = t0;
        for (int i = 0; i < dim; i ++) {
            xrms[i] = 0.0;
            overallxrms[i] = 0.0;
        }
        // Check the mode of the bifurcation
        if (bifmode == 1) {
            // Reset Initial conditions in each bifurcation step
            for (int i = 0; i < dim; i++) {
                x[i] = IC[i];
            }
        }
        // Reset initial values to xmax and xmin based on initial values of x
        for (int i = 0; i < dim; i++) {
            xmax[i] = x[i];
            xmin[i] = x[i];
        }
        // Vary timestep if varpar = par[0]
        h = (2 * pi) / (ndiv * par[0]);         // par[0] = OMEGA
        // Call Runge-Kutta 4th order integrator n = np * ndiv times
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < ndiv; j++) {
                rk4(dim, x, t, h, par, f, edosys);
                t = t + h;
                // Apply poincare map at permanent regime
                if (i > trans) {
                    // Get max and min values at permanent regime
                    for (int q = 0; q < dim; q++) {
                        max_value(x[q], &xmax[q]);
                        min_value(x[q], &xmin[q]);
                    }
                    // Choose any point in the trajectory for plane placement
                    if (j == 1) {
                        // Print the result in output file
                        write_results(output_poinc_file, dim, par[parindex], x, xmin, xmax, nrms, rmsindex, xrms, overallxrms, 1);
                    }
                    // Accumulate squared values to RMS computation in permanent regime
                    for (int q = 0; q < nrms; q++) {
                        xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], x[rmsindex[q]], N, 0);
                    }
                }
                // Accumulate squared values to RMS computation for all time domain
                for (int q = 0; q < nrms; q++) {
                    overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], x[rmsindex[q]], N, 0);
                }
            }
        }
        // Compute RMS values of state variables
        for (int q = 0; q < nrms; q++) {
            xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], x[rmsindex[q]], N, 1);
            overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], x[rmsindex[q]], N, 1);
        }
        // Print results in output file
        write_results(output_file, dim, par[parindex], x, xmin, xmax, nrms, rmsindex, xrms, overallxrms, 3);
        // Progress Monitor
        if (parrange[2] > 100) {
            if (k % 50 == 0) {
                progress_bar(0, par[parindex], parrange[0], parrange[1]);
            }
            if (k == parrange[2] - 1) {
                progress_bar(9, par[parindex], parrange[0], parrange[1]);
            }
        } else {
            progress_bar(0, par[parindex], parrange[0], parrange[1]);
        }
        
    }
    // Free Memory
    free(f); free(IC); free(xmax); free(xmin); free(xrms); free(overallxrms);
}

void EH_full_bifurcation_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int maxper, double t,
                                  double **x, int parindex, double *parrange, double *par, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *),
                                  void (*write_results)(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *LE, int attractor, double **poinc, int diffattrac, int nrms, int *rmsindex, double *xrms, double *overallxrms, int mode), int bifmode) {
    // Allocate memory for x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Allocate memory for vectors necessary for lyapunov exponents calculation
    double *cum = malloc(dim * sizeof *cum);            // Cumulative Vector
    double *lambda = malloc(dim *sizeof *lambda);       // Lyapunov Exponents vector
    double *s_cum = malloc(dim * sizeof *s_cum);        // Short Cumulative Vector
    double *s_lambda = malloc(dim * sizeof *s_lambda);  // Short Lyapunov Exponents Vector
    double *znorm = malloc(dim * sizeof *znorm);        // Norm of Vectors
    double *gsc = malloc((dim - 1) * sizeof *gsc);      // Inner Products Vector
    // Store Initial Conditions
    double t0 = t;
    double *IC = malloc(dim * sizeof *IC);
    for (int i = 0; i < dim; i++) {
        IC[i] = (*x)[i];
    }
    // Declare memory to store min and max values
    double *xmax = malloc(dim * sizeof *xmax);
    double *xmin = malloc(dim * sizeof *xmin);
    // Allocate RMS variables
    double *xrms = malloc(dim * sizeof *xrms);
    double *overallxrms = malloc(dim * sizeof *overallxrms);
    // Mumber of integration steps
    int N = np*ndiv;
    // Declare timestep, final time, short initial time and pi
    double h, tf, s_T0;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control parameters
    double varstep;        
    varstep = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    // Numerical control parameters
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations 
    // Declare vector and allocate memory to store poincare map values: poinc[number of permanent regime forcing periods][dimension original system]
    double **poinc = malloc((np - trans) * sizeof **poinc);
    for (int i = 0; i < np - trans; i++) {
        poinc[i] = malloc(dim * sizeof **poinc);
    }
    // Declare vector to store the chosen Lyapunov Exponents to determine the attractor
    double *LE = malloc(dim * sizeof *LE);
    // Declare variable to store attractor
    int attrac = 0;
    // Declare vector for temporary storage of periodicity values to check if all directions are equal
    int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
    // Declare variable to store if all directions have the same periodicity (0 if all directions are equal, 1 if all directions arent equal)
    int diffAttrac = -1;
    // Prepare x vector to include perturbed values
    realloc_vector(x, ndim);
    // Make the header of output files
    write_results(output_poinc_file, dim, np, trans, par[parindex], (*x), xmin, xmax, LE, attrac, poinc, diffAttrac, nrms, rmsindex, xrms, overallxrms, 0);
    write_results(output_file, dim, np, trans, par[parindex], (*x), xmin, xmax, LE, attrac, poinc, diffAttrac, nrms, rmsindex, xrms, overallxrms, 2);
    // Starts to increment bifurcation control parameter
    for (int k = 0; k < (int)parrange[2]; k++) {
        par[parindex] = parrange[0] + k*varstep; // Increment value
        // Reset Variables
        t = t0;
        for (int i = 0; i < dim; i++) {
            lambda[i] = 0.0;
            s_lambda[i] = 0.0;
            LE[i] = 0.0;
            xrms[i] = 0.0;
            overallxrms[i] = 0.0;
        }
        // Check the mode of the bifurcation
        if (bifmode == 1) {
            // Reset Initial conditions in each bifurcation step
            for (int i = 0; i < dim; i++) {
                (*x)[i] = IC[i];
            }
        }
        // Reset initial values to xmax and xmin based on initial values of x
        for (int i = 0; i < dim; i++) {
            xmax[i] = (*x)[i];
            xmin[i] = (*x)[i];
        }
        // Vary timestep if varpar = par[0], varying also final time and short initial time
        h = (2 * pi) / (ndiv * par[0]);              // par[0] = OMEGA
        tf = h*np*ndiv;                              // Final time
        s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time
        // Assign initial perturbation
        perturb_wolf(x, dim, ndim, &cum, &s_cum);
        // Call Runge-Kutta 4th order integrator N = nP * nDiv times
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < ndiv; j++) {
                rk4(ndim, *x, t, h, par, f, edosys);
                lyapunov_wolf(x, t, h, dim, ndim, s_T0, &cum, &s_cum, &lambda, &s_lambda, &znorm, &gsc);
                t = t + h;
                // Apply poincare map at permanent regime
                if (i >= trans) {
                    // Get max and min values at permanent regime
                    for (int q = 0; q < dim; q++) {
                        max_value((*x)[q], &xmax[q]);
                        min_value((*x)[q], &xmin[q]);
                    }
                    // Accumulate squared values to RMS computation in permanent regime
                    for (int q = 0; q < nrms; q++) {
                        xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], (*x)[rmsindex[q]], N, 0);
                    }
                    // Choose any point in the trajectory for poincare section placement
                    if (j == 1) {
                        // Stores poincare values in poinc[np - trans][dim] vector
                        for (int p = 0; p < dim; p++) {
                            poinc[i - trans][p] = (*x)[p];
                        }
                    }
                }
                // Accumulate squared values to RMS computation for all time domain
                for (int q = 0; q < nrms; q++) {
                    overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], (*x)[rmsindex[q]], N, 0);
                }
            }
        }
        // Compute RMS values of state variables
        for (int q = 0; q < nrms; q++) {
            xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
            overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
        }
        // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
        store_LE(dim, lambda, s_lambda, LE);
        // Verify the type of motion of the system
        attrac = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
        // Write results in file
        write_results(output_poinc_file, dim, np, trans, par[parindex], (*x), xmin, xmax, LE, attrac, poinc, diffAttrac, nrms, rmsindex, xrms, overallxrms, 1);
        write_results(output_file, dim, np, trans, par[parindex], (*x), xmin, xmax, LE, attrac, poinc, diffAttrac, nrms, rmsindex, xrms, overallxrms, 3);
        // Progress Monitor
        if (parrange[2] > 100) {
            if (k % 50 == 0) {
                progress_bar(0, par[parindex], parrange[0], parrange[1]);
            }
            if (k == parrange[2] - 1) {
                progress_bar(0, par[parindex], parrange[0], parrange[1]);
            }
        } else {
            progress_bar(0, par[parindex], parrange[0], parrange[1]);
        }
    }
    // Free memory    
    free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
    free(znorm); free(gsc); free(LE); free(tmp_attrac); free(IC);
    free(xmax); free(xmin);
    for (int i = 0; i < np - trans; i++) {
        free(poinc[i]);
    }
    free(poinc);
    free(xrms); free(overallxrms);
}