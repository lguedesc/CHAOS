#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "iofiles.h"
#include "nldyn.h"
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1 
#endif

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

// Solutions
void OS_rk4_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, int nrms, int *rmsindex, double **xrms, double **overallxrms, 
                     double **xmin, double **xmax, double **overallxmin, double **overallxmax, void (*edosys)(int, double *, double, double *, double *),  
                     int ncustomvalues, char ***customnames, double **customvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex,
                     void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode)) {
    // Maximum length of custom names, if there is custom calculations
    size_t nchars = 20;
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
    // Declare memory to store min, max, overall min and overall max values
    (*xmax) = malloc(dim * sizeof **xmax);
    (*xmin) = malloc(dim * sizeof **xmin);
    (*overallxmax) = malloc(dim * sizeof **xmin);
    (*overallxmin) = malloc(dim * sizeof **xmin);
    // Initialize min, max, overall min and overall max vectors
    for (int i = 0; i < dim; i ++) {
        (*xmin)[i] = 0.0;
        (*overallxmin)[i] = x[i];
        (*xmax)[i] = 0.0;
        (*overallxmax)[i] = x[i];
    }
    if (ncustomvalues > 0) {
        // Allocate variable to store custom values
        (*customvalues) = malloc(ncustomvalues * sizeof **customvalues);
        // Initialize customvalues[ncustomvalues]
        for (int i = 0; i < ncustomvalues; i++) {
            (*customvalues)[i] = 0.0;
        }
        // Allocate variable to store the names of the custom values
        (*customnames) = malloc(ncustomvalues * sizeof *customnames);
        for (int i = 0; i < ncustomvalues; i++) {
            (*customnames)[i] = malloc(nchars * sizeof *customnames);
        }
    } 
    // Mumber of integration steps
    int N = np*ndiv;
    // Percentage of integration steps considered steady state regime 
    double steadystateperc = 1 - ((double)trans/(double)np);
    // Check if there is any custom names to be inserted on the header of the output file
    if (ncustomvalues > 0) {
        customfunc(x, par, t, (*xrms), (*xmin), (*xmax), N, steadystateperc, ncustomvalues, (*customnames), nchars, (*customvalues), 0);
    }
    // Make the header of the output file
    EH_write_timeseries_results(output_file, dim, t, x, ncustomvalues, (*customnames), (*customvalues), nprintf, printfindex, 1);
    // Call Runge-Kutta 4th order integrator n = np * ndiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(dim, x, t, h, par, f, edosys);
            t = t + h;
            // Steady State Regime Calculations
            if (i >= trans) {
                // Get max and min values at steady state regime
                for (int q = 0; q < dim; q++) {
                    // Initialize xmax[dim] and xmin[dim] with first values of x[dim] at steady state regime
                    if (i == trans && j == 0) {
                        (*xmax)[q] = x[q];
                        (*xmin)[q] = x[q];
                    }
                    max_value(x[q], &(*xmax)[q]);
                    min_value(x[q], &(*xmin)[q]);
                }
                // Accumulate squared values to RMS computation in steady state regime (if there is RMS calculations to be performed)
                if (nrms > 0) {
                    for (int q = 0; q < nrms; q++) {
                        (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], x[rmsindex[q]], N, 0);
                    }
                }
                // Perform custom calculations in steady state regime (if there is calculations to be done)
                if (ncustomvalues > 0) {
                    customfunc(x, par, t, (*xrms), (*xmin), (*xmax), N, steadystateperc, ncustomvalues, (*customnames), nchars, (*customvalues), 1);    
                }
            }
            // Get overall max and min values
            for (int q = 0; q < dim; q++) {
                    //printf("overallxmax[%d] = %lf\n", q, overallxmax[q]);
                    max_value(x[q], &(*overallxmax)[q]);
                    min_value(x[q], &(*overallxmin)[q]);
            }
            // Accumulate squared values to RMS computation for all time domain (if there is RMS calculations to be performed)
            if (nrms > 0) {
                for (int q = 0; q < nrms; q++) {
                    (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], x[rmsindex[q]], N, 0);
                }
            }
            // Perform custom calculations over the entire time series (transient + steady state), if there is calculations to be done
            if (ncustomvalues > 0) {
                customfunc(x, par, t, (*xrms), (*xmin), (*xmax), N, steadystateperc, ncustomvalues, (*customnames), nchars, (*customvalues), 2);
            }
            // Write results in output file
            EH_write_timeseries_results(output_file, dim, t, x, ncustomvalues, (*customnames), (*customvalues), nprintf, printfindex, 2);
        }
    }
    // Compute RMS values of state variables (if there is RMS calculations to be performed)
    if (nrms > 0) {
        for (int q = 0; q < nrms; q++) {
            (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], x[rmsindex[q]], N, 1);
            (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], x[rmsindex[q]], N, 1);
        }
    }
    // Perform custom calculations at the end of the time series (if there is calculations to be done)
    if (ncustomvalues > 0) {
        customfunc(x, par, t, (*xrms), (*xmin), (*xmax), N, steadystateperc, ncustomvalues, (*customnames), nchars, (*customvalues), 3);
    }
    // Free Memory
    free(f); 
}

void OS_full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, 
                                 int nrms, int *rmsindex, double **xrms, double **overallxrms, double **xmin, double **xmax, double **overallxmin, double **overallxmax,
                                 void (*edosys)(int, double *, double, double *, double *), 
                                 int ncustomvalues, char ***customnames, double **customvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex,
                                 void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode)) {
    // Maximum length of custom names, if there is custom calculations
    size_t nchars = 20;
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
    // Declare memory to store min, max, overall min and overall max values
    (*xmax) = malloc(dim * sizeof **xmax);
    (*xmin) = malloc(dim * sizeof **xmin);
    (*overallxmax) = malloc(dim * sizeof **overallxmax);
    (*overallxmin) = malloc(dim * sizeof **overallxmin);
    // Initialize min, max, overall min and overall max vectors
    for (int i = 0; i < dim; i ++) {
        (*xmin)[i] = 0.0;
        (*overallxmin)[i] = (*x)[i];
        (*xmax)[i] = 0.0;
        (*overallxmax)[i] = (*x)[i];
    }
    if (ncustomvalues > 0) {
        // Allocate variable to store custom values
        (*customvalues) = malloc(ncustomvalues * sizeof **customvalues);
        // Initialize customvalues[ncustomvalues]
        for (int i = 0; i < ncustomvalues; i++) {
            (*customvalues)[i] = 0.0;
        }
        // Allocate variable to store the names of the custom values
        (*customnames) = malloc(ncustomvalues * sizeof *customnames);
        for (int i = 0; i < ncustomvalues; i++) {
            (*customnames)[i] = malloc(nchars * sizeof *customnames);
        }
    }
    // Mumber of integration steps
    int N = np*ndiv;
    // Percentage of integration steps considered steady state regime 
    double steadystateperc = 1 - ((double)trans/(double)np);
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
    // Check if there is any custom names to be inserted on the header of the output file
    if (ncustomvalues > 0) {
        customfunc((*x), par, t, (*xrms), (*xmin), (*xmax), N, steadystateperc, ncustomvalues, (*customnames), nchars, (*customvalues), 0);
    }
    // Make the header of output files
    EH_write_ftimeseries_results(output_ftimeseries_file, dim, t, (*x), lambda, s_lambda, ncustomvalues, (*customnames), (*customvalues), nprintf, printfindex, 1);
    write_results(output_poinc_file, dim, t, (*x), 0);
    // Call Runge-Kutta 4th order integrator N = nP * nDiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(ndim, *x, t, h, par, f, edosys);
            lyapunov_wolf(x, t, h, dim, ndim, s_T0, &cum, &s_cum, &lambda, &s_lambda, &znorm, &gsc);
            t = t + h;
            // Steady State Regime Calculations
            if (i >= trans) {
                // Choose any point in the trajectory for poincare section placement
                if (j == 1) {
                    // Stores poincare values in poinc[np - trans][dim] vector
                    for (int p = 0; p < dim; p++) {
                        poinc[i - trans][p] = (*x)[p];
                    }
                    write_results(output_poinc_file, dim, t, (*x), 2);
                }
                // Get max and min values at permanent regime
                for (int q = 0; q < dim; q++) {
                    // Initialize xmax[dim] and xmin[dim] with first values of x[dim] at permanent regime
                    if (i == trans && j == 0) {
                        (*xmax)[q] = (*x)[q];
                        (*xmin)[q] = (*x)[q];
                    }
                    max_value((*x)[q], &(*xmax)[q]);
                    min_value((*x)[q], &(*xmin)[q]);
                }
                // Accumulate squared values to RMS computation in permanent regime
                if (nrms > 0) {
                    for (int q = 0; q < nrms; q++) {
                        (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 0);
                    }
                }
                // Perform custom calculations in steady state regime (if there is calculations to be done)
                if (ncustomvalues > 0) {
                    customfunc((*x), par, t, (*xrms), (*xmin), (*xmax), N, steadystateperc, ncustomvalues, (*customnames), nchars, (*customvalues), 1);    
                }
            }
            // Get overall max and min values
            for (int q = 0; q < dim; q++) {
                    //printf("overallxmax[%d] = %lf\n", q, overallxmax[q]);
                    max_value((*x)[q], &(*overallxmax)[q]);
                    min_value((*x)[q], &(*overallxmin)[q]);
            }
            // Accumulate squared values to RMS computation for all time domain
            if (nrms > 0) {
                for (int q = 0; q < nrms; q++) {
                    (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 0);
                }
            }
            // Perform custom calculations over the entire time series (transient + steady state), if there is calculations to be done
            if (ncustomvalues > 0) {
                customfunc((*x), par, t, (*xrms), (*xmin), (*xmax), N, steadystateperc, ncustomvalues, (*customnames), nchars, (*customvalues), 2);
            }       
            // Write Results in output file
            EH_write_ftimeseries_results(output_ftimeseries_file, dim, t, (*x), lambda, s_lambda, ncustomvalues, (*customnames), (*customvalues), nprintf, printfindex, 2);
        }
    }
    // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
    store_LE(dim, lambda, s_lambda, LE);
    // Declare vector for temporary storage of periodicity values to check if all directions are equal
    int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
    // Declare variable to flag if all directions present same periodicity or not (0 = all the same, 1 = not the same)
    int diffAttrac = -1;
    // Verify the type of motion of the system
    (*attrac) = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
    // Compute RMS values of state variables
    if (nrms > 0) {
        for (int q = 0; q < nrms; q++) {
            (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
            (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
        }
    }
    // Perform "end" type custom calculations if there is calculations to be done
    if (ncustomvalues > 0) {
        customfunc((*x), par, t, (*xrms), (*xmin), (*xmax), N, steadystateperc, ncustomvalues, (*customnames), nchars, (*customvalues), 3);
    }
    // Free memory    
    free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
    free(znorm); free(gsc); free(LE); free(tmp_attrac);
    for (int i = 0; i < np - trans; i++) {
        free(poinc[i]);
    }
    free(poinc);
}

void OS_bifurc_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, double t, double *x, int parindex, 
                        double *parrange, double *par, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                        int ncustomvalues, int nprintf, int *printfindex,
                        void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                        int mode), int bifmode) {
    // Maximum length of custom names, if there is custom calculations
    size_t nchars = 20;
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
    double *overallxmin = malloc(dim * sizeof *overallxmin);
    double *overallxmax = malloc(dim * sizeof *overallxmax);
    // Allocate RMS variables
    double *xrms = malloc(dim * sizeof *xrms);
    double *overallxrms = malloc(dim * sizeof *overallxrms);
    // Allocate variable to store custom values
    double *customvalues = malloc(ncustomvalues * sizeof *customvalues);
    // Allocate variable to store the names of the custom values
    char **customnames = malloc(ncustomvalues * sizeof *customnames);
    for (int i = 0; i < ncustomvalues; i++) {
        customnames[i] = malloc(nchars * sizeof *customnames);
    }
    // Mumber of integration steps
    int N = np*ndiv;
    // Percentage of integration steps considered steady state regime 
    double steadystateperc = 1 - ((double)trans/(double)np);
    // Declare timestep and pi
    double h;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control parameters
    double varstep;        
    varstep = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    // Check if there is any custom names to be inserted on the header of the output file
    if (ncustomvalues > 0) {
        customfunc(x, par, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 0);
    }
    // Make the header of the output file
    EH_write_bifurc_results(output_poinc_file, dim, par[parindex], x, xmin, xmax, overallxmin, overallxmax, nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, 0);    
    EH_write_bifurc_results(output_file, dim, par[parindex], x, xmin, xmax, overallxmin, overallxmax, nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, 2);
    // Starts sweep the control parameter
    for (int k = 0; k < (int)parrange[2]; k++) {
        par[parindex] = parrange[0] + k*varstep; // Increment value
        // Check the mode of the bifurcation
        if (bifmode == 1) {
            // Reset Initial conditions in each bifurcation step
            for (int i = 0; i < dim; i++) {
                x[i] = IC[i];
            }
        }
        // Reset Variables
        t = t0;
        for (int i = 0; i < dim; i ++) {
            xmin[i] = 0.0;
            overallxmin[i] = x[i];
            xmax[i] = 0.0;
            overallxmax[i] = x[i];
            xrms[i] = 0.0;
            overallxrms[i] = 0.0;
        }
        if (ncustomvalues > 0) {
            for (int i = 0; i < ncustomvalues; i++) {
                customvalues[i] = 0.0;
            }
        }
        // Vary timestep if varpar = par[0]
        h = (2 * pi) / (ndiv * par[0]);         // par[0] = OMEGA
        // Call Runge-Kutta 4th order integrator n = np * ndiv times
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < ndiv; j++) {
                rk4(dim, x, t, h, par, f, edosys);
                t = t + h;
                // Apply poincare map at permanent regime
                if (i >= trans) {
                    // Get max and min values at permanent regime
                    for (int q = 0; q < dim; q++) {
                        // Initialize xmax[dim] and xmin[dim] with first values of x[dim] at permanent regime
                        if (i == trans && j == 0) {
                            xmax[q] = x[q];
                            xmin[q] = x[q];
                        }
                        max_value(x[q], &xmax[q]);
                        min_value(x[q], &xmin[q]);
                    }
                    // Choose any point in the trajectory for plane placement
                    if (j == 1) {
                        // Print the result in output file
                        EH_write_bifurc_results(output_poinc_file, dim, par[parindex], x, xmin, xmax, overallxmin, overallxmax, nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, 1);
                    }
                    // Accumulate squared values to RMS computation in permanent regime
                    if (nrms > 0) {
                        for (int q = 0; q < nrms; q++) {
                            xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], x[rmsindex[q]], N, 0);
                        }
                    }
                    // Perform custom calculations in steady state regime (if there is calculations to be done)
                    if (ncustomvalues > 0) {
                        customfunc(x, par, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 1);    
                    }
                }
                // Get overall max and min values
                for (int q = 0; q < dim; q++) {
                    max_value(x[q], &overallxmax[q]);
                    min_value(x[q], &overallxmin[q]);
                }
                // Accumulate squared values to RMS computation for all time domain
                if (nrms > 0) {
                    for (int q = 0; q < nrms; q++) {
                        overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], x[rmsindex[q]], N, 0);
                    }
                }
                // Perform custom calculations over the entire time series (transient + steady state), if there is calculations to be done
                if (ncustomvalues > 0) {
                    customfunc(x, par, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 2);
                }
            }
        }
        // Compute RMS values of state variables
        if (nrms > 0) {
            for (int q = 0; q < nrms; q++) {
                xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], x[rmsindex[q]], N, 1);
                overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], x[rmsindex[q]], N, 1);
            }
        }
        // Perform custom calculations at the end of the time series (if there is calculations to be done)
        if (ncustomvalues > 0) {
            customfunc(x, par, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 3);
        }
        // Print results in output file
        EH_write_bifurc_results(output_file, dim, par[parindex], x, xmin, xmax, overallxmin, overallxmax, nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, 3);
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
    free(overallxmin); free(overallxmax); free(customvalues);
    for (int i = 0; i < ncustomvalues; i++) {
            free(customnames[i]);
        }
    free(customnames);
}

void OS_full_bifurcation_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int maxper, double t,
                                  double **x, int parindex, double *parrange, double *par, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                                  int ncustomvalues, int nprintf, int *printfindex,
                                  void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                                  int mode), int bifmode) {
    // Maximum length of custom names, if there is custom calculations
    size_t nchars = 20;
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
    double *overallxmin = malloc(dim * sizeof *overallxmin);
    double *overallxmax = malloc(dim * sizeof *overallxmax);
    // Allocate RMS variables
    double *xrms = malloc(dim * sizeof *xrms);
    double *overallxrms = malloc(dim * sizeof *overallxrms);
    // Allocate variable to store custom values
    double *customvalues = malloc(ncustomvalues * sizeof *customvalues);
    // Allocate variable to store the names of the custom values
    char **customnames = malloc(ncustomvalues * sizeof *customnames);
    for (int i = 0; i < ncustomvalues; i++) {
        customnames[i] = malloc(nchars * sizeof *customnames);
    }
    // Mumber of integration steps
    int N = np*ndiv;
    // Percentage of integration steps considered steady state regime 
    double steadystateperc = 1 - ((double)trans/(double)np);
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
    // Check if there is any custom names to be inserted on the header of the output file
    if (ncustomvalues > 0) {
        customfunc((*x), par, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 0);
    }
    // Make the header of output files
    EH_write_fbifurc_results(output_poinc_file, dim, np, trans, par[parindex], (*x), xmin, xmax, overallxmin, overallxmax, LE, attrac, poinc, diffAttrac,
                             nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, 0);
    EH_write_fbifurc_results(output_file, dim, np, trans, par[parindex], (*x), xmin, xmax, overallxmin, overallxmax, LE, attrac, poinc, diffAttrac,
                             nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, 2);
    // Starts to increment bifurcation control parameter
    for (int k = 0; k < (int)parrange[2]; k++) {
        par[parindex] = parrange[0] + k*varstep; // Increment value
        // Check the mode of the bifurcation
        if (bifmode == 1) {
            // Reset Initial conditions in each bifurcation step
            for (int i = 0; i < dim; i++) {
                (*x)[i] = IC[i];
            }
        }
        // Reset Variables
        t = t0;
        for (int i = 0; i < dim; i++) {
            lambda[i] = 0.0;
            s_lambda[i] = 0.0;
            LE[i] = 0.0;
            xmin[i] = 0.0;
            overallxmin[i] = (*x)[i];
            xmax[i] = 0.0;
            overallxmax[i] = (*x)[i];
            xrms[i] = 0.0;
            overallxrms[i] = 0.0;
        }
        if (ncustomvalues > 0) {
            for (int i = 0; i < ncustomvalues; i++) {
                customvalues[i] = 0.0;
            }
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
                        // Initialize xmax[dim] and xmin[dim] with first values of x[dim] at permanent regime
                        if (i == trans && j == 0) {
                            xmax[q] = (*x)[q];
                            xmin[q] = (*x)[q];
                        }
                        max_value((*x)[q], &xmax[q]);
                        min_value((*x)[q], &xmin[q]);
                    }
                    // Accumulate squared values to RMS computation in permanent regime
                    if (nrms > 0) {
                        for (int q = 0; q < nrms; q++) {
                            xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], (*x)[rmsindex[q]], N, 0);
                        }
                    }
                    // Perform custom calculations in steady state regime (if there is calculations to be done)
                    if (ncustomvalues > 0) {
                        customfunc((*x), par, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 1);    
                    }
                    // Choose any point in the trajectory for poincare section placement
                    if (j == 1) {
                        // Stores poincare values in poinc[np - trans][dim] vector
                        for (int p = 0; p < dim; p++) {
                            poinc[i - trans][p] = (*x)[p];
                        }
                    }
                }
                // Get overall max and min values
                for (int q = 0; q < dim; q++) {
                    max_value((*x)[q], &overallxmax[q]);
                    min_value((*x)[q], &overallxmin[q]);
                }
                // Accumulate squared values to RMS computation for all time domain
                if (nrms > 0) {
                    for (int q = 0; q < nrms; q++) {
                        overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], (*x)[rmsindex[q]], N, 0);
                    }
                }
                // Perform custom calculations over the entire time series (transient + steady state), if there is calculations to be done
                if (ncustomvalues > 0) {
                    customfunc((*x), par, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 2);
                }
            }
        }
        // Compute RMS values of state variables
        for (int q = 0; q < nrms; q++) {
            xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
            overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
        }
        // Perform custom calculations at the end of the time series (if there is calculations to be done)
        if (ncustomvalues > 0) {
            customfunc((*x), par, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 3);
        }
        // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
        store_LE(dim, lambda, s_lambda, LE);
        // Verify the type of motion of the system
        attrac = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
        // Write results in file
        EH_write_fbifurc_results(output_poinc_file, dim, np, trans, par[parindex], (*x), xmin, xmax, overallxmin, overallxmax, LE, attrac, poinc, diffAttrac,
                             nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, 1);
        EH_write_fbifurc_results(output_file, dim, np, trans, par[parindex], (*x), xmin, xmax, overallxmin, overallxmax, LE, attrac, poinc, diffAttrac,
                             nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, 3);
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
    free(xmax); free(xmin); free(overallxmax); free(overallxmin);
    for (int i = 0; i < np - trans; i++) {
        free(poinc[i]);
    }
    free(poinc);
    free(xrms); free(overallxrms);
    free(customvalues);
    for (int i = 0; i < ncustomvalues; i++) {
            free(customnames[i]);
        }
    free(customnames);
}

void OS_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                    int indexX, int indexY, double *parrange, double *par, int npar, int nrms, int *rmsindex,
                                    void (*edosys)(int, double *, double, double *, double *),
                                    int ncustomvalues, int nprintf, int *printfindex,
                                    void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                                    int mode), int bifmode) {
    // Maximum length of custom names, if there is custom calculations
    size_t nchars = 20;
    // Declare matrix do store results
    int pixels = parrange[2]*parrange[5];  // Number of results
    double **results = malloc(pixels * sizeof **results);
    for (int i = 0; i < pixels; i++) {
        //results[i] = malloc((4 + (3*dim)) * sizeof **results);  // 4 for params, attrac, diffratac/ dim for xmax/ dim for xmin / dim for LE
        results[i] = malloc((4 + (4*dim) + (2*nrms) + nprintf) * sizeof **results);
    }
    // Declare rk4 timestep, final time, short initial time, pi and varstep
    double h;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control parameters
    double varstep[2];        
    varstep[0] = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    varstep[1] = (parrange[4] - parrange[3])/(parrange[5] - 1); // -1 in the denominator ensures the input resolution
    // Declare variable to flag if all directions present same periodicity or not (0 = all the same, 1 = not the same)
    int diffAttrac = -1;
    // Declare variable to store attractor
    int attrac;
    // Allocate variable to store the names of the custom values
    char **customnames = malloc(ncustomvalues * sizeof *customnames);
    for (int i = 0; i < ncustomvalues; i++) {
        customnames[i] = malloc(nchars * sizeof *customnames);
    }
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, bifmode, ndiv, np, trans, maxper, varstep, npar, results, edosys, customfunc, customnames, nchars, pi, rmsindex, nrms, printfindex, nprintf, ncustomvalues) \
                                       private(t, h, attrac) \
                                       firstprivate(x, indexX, indexY, parrange, par, diffAttrac)
    {   
        //Get number of threads
        int ID = omp_get_thread_num();
        // Allocate memory for x` = f(x)
        double *f = malloc(dim * sizeof *f);
        // Declare vector and allocate memory to store poincare map values: poinc[number of permanent regime forcing periods][dimension original system]
        double **poinc = malloc((np - trans) * sizeof **poinc);
        for (int i = 0; i < np - trans; i++) {
            poinc[i] = malloc(dim * sizeof **poinc);
        }
        // Declare vector for temporary storage of periodicity values to check if all directions are equal
        int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
        // Store Initial Conditions
        double t0 = t;
        double *IC = malloc(dim * sizeof *IC);
        for (int i = 0; i < dim; i++) {
            IC[i] = (*x)[i];
        }
        // Declare memory to store min and max values
        double *xmax = malloc(dim * sizeof *xmax);
        double *xmin = malloc(dim * sizeof *xmin);
        double *overallxmin = malloc(dim * sizeof *overallxmin);
        double *overallxmax = malloc(dim * sizeof *overallxmax);
        // Allocate RMS variables
        double *xrms = malloc(dim * sizeof *xrms);
        double *overallxrms = malloc(dim * sizeof *overallxrms);
        // Allocate variable to store custom values
        double *customvalues = malloc(ncustomvalues * sizeof *customvalues);
        // Mumber of integration steps
        int N = np*ndiv;
        // Percentage of integration steps considered steady state regime 
        double steadystateperc = 1 - ((double)trans/(double)np);
        // Convert function arguments as local (private) variables
        double *X = convert_argument_to_private(*x, dim);
        double *PAR = convert_argument_to_private(par, npar);
        // Check if there is any custom names to be inserted on the header of the output file
        if (ncustomvalues > 0) {
            customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 0);
        }
        // Index to identify position to write results
        int index;                                          
        // Declare counter for parallelized loop
        int k, m;
        #pragma omp for schedule(static) private(k, m)
        // Starts the parallel loop for Y control parameter
        for (k = 0; k < (int)parrange[5]; k++) {
            PAR[indexY] = parrange[3] + k*varstep[1]; // Increment value
            // Reset Initial conditions for the beggining of a horizontal line
            for (int i = 0; i < dim; i++) {
                X[i] = IC[i];
            }
            // Starts the loop for X control parameter
            for (m = 0; m < (int)parrange[2]; m++) {
                PAR[indexX] = parrange[0] + m*varstep[0]; // Increment Value
                // Check the mode of the diagram
                if (bifmode == 1) {
                    // Reset Initial conditions in each diagram step
                    for (int i = 0; i < dim; i++) {
                        X[i] = IC[i];
                    }
                }
                // Reset Variables
                t = t0;
                for (int i = 0; i < dim; i++) {
                    xrms[i] = 0.0;
                    overallxrms[i] = 0.0;
                    xmin[i] = 0.0;
                    overallxmin[i] = X[i];
                    xmax[i] = 0.0;
                    overallxmax[i] = X[i];
                }
                if (ncustomvalues > 0) {
                    for (int i = 0; i < ncustomvalues; i++) {
                        customvalues[i] = 0.0;
                    }
                }
                // Vary timestep if varpar = par[0], varying also final time and short initial time
                h = (2 * pi) / (ndiv * PAR[0]);              // par[0] = OMEGA
                // Call Runge-Kutta 4th order integrator N = nP * nDiv times
                for (int i = 0; i < np; i++) {
                    for (int j = 0; j < ndiv; j++) {
                        rk4(dim, X, t, h, PAR, f, edosys);
                        t = t + h;
                        // Apply poincare map at permanent regime
                        if (i >= trans) {
                            // Get max and min values at permanent regime
                            for (int q = 0; q < dim; q++) {
                                // Initialize xmax[dim] and xmin[dim] with first values of x[dim] at permanent regime
                                if (i == trans && j == 0) {
                                    xmax[q] = X[q];
                                    xmin[q] = X[q];
                                }
                                max_value(X[q], &xmax[q]);
                                min_value(X[q], &xmin[q]);
                            }
                            // Accumulate squared values to RMS computation in permanent regime
                            if (nrms > 0) {
                                for (int q = 0; q < nrms; q++) {
                                    xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], X[rmsindex[q]], N, 0);
                                }
                            }
                            // Perform custom calculations in steady state regime (if there is calculations to be done)
                            if (ncustomvalues > 0) {
                                customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 1);
                            }
                            // Choose any point in the trajectory for poincare section placement
                            if (j == 1) {
                                // Stores poincare values in poinc[np - trans][dim] vector
                                for (int p = 0; p < dim; p++) {
                                    poinc[i - trans][p] = X[p];
                                }
                            }
                        }
                        // Get overall max and min values
                        for (int q = 0; q < dim; q++) {
                            max_value(X[q], &overallxmax[q]);
                            min_value(X[q], &overallxmin[q]);
                        }
                        // Accumulate squared values for RMS computation for all time domain
                        if (nrms > 0) {
                            for (int q = 0; q < nrms; q++) {
                                overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], X[rmsindex[q]], N, 0);
                            }
                        }
                        // Perform custom calculations over the entire time series (transient + steady state), if there is calculations to be done
                        if (ncustomvalues > 0) {
                            customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 2);
                        }
                    }
                }
                // Compute RMS values of state variables
                if (nrms > 0) {
                    for (int q = 0; q < nrms; q++) {
                        xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], X[rmsindex[q]], N, 1);
                        overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], X[rmsindex[q]], N, 1);
                    }
                }
                // Perform custom calculations at the end of the time series (if there is calculations to be done)
                if (ncustomvalues > 0) {
                    customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 3);
                }
                // Verify the type of motion of the system
                attrac = check_periodicity(dim, np, poinc, trans, tmp_attrac, &diffAttrac, maxper);
                // Write results in matrix
                index = (int)parrange[2]*k + m;
                results[index][0] = PAR[indexY];
                results[index][1] = PAR[indexX];
                results[index][2] = (double)attrac;
                results[index][3] = (double)diffAttrac;
                for (int r = 4; r < dim + 4; r++) {
                    results[index][r] = xmax[r-4];
                }
                for (int r = dim + 4; r < 4 + (2*dim); r++) {
                    results[index][r] = xmin[r - 4 - dim];
                }
                for (int r = 4 + (2*dim); r < 4 + (3*dim); r++) {
                    results[index][r] = overallxmax[r - 4 - (2*dim)];
                }
                for (int r = 4 + (3*dim); r < 4 + (4*dim); r++) {
                    results[index][r] = overallxmin[r - 4 - (3*dim)];
                }
                for (int r = 4 + (4*dim); r < 4 + (4*dim) + nrms; r++) {
                    results[index][r] = xrms[rmsindex[r - 4 - (4*dim)]];
                }
                for (int r = 4 + (4*dim) + nrms; r < 4 + (4*dim) + (2*nrms); r++) {
                    results[index][r] = overallxrms[rmsindex[r - 4 - (4*dim) - nrms]];
                }
                for (int r = 4 + (4*dim) + (2*nrms); r < 4 + (4*dim) + (2*nrms) + nprintf; r ++) {
                    results[index][r] = customvalues[printfindex[r - 4 - (4*dim) - (2*nrms)]];
                }
            }
            // Progress Monitor
            if (ID == 0) {
                progress_bar(0, PAR[indexY], parrange[3], (parrange[4] - varstep[1])/omp_get_num_threads());
                if (k == ((int)parrange[5] - 1)/omp_get_num_threads() ) {
                    progress_bar(1, PAR[indexY], parrange[3], (parrange[4] - varstep[1])/omp_get_num_threads());
                }
            }
        }
        // Free memory    
        free(f);
        free(tmp_attrac); free(IC);
        free(xmax); free(xmin); free(overallxmax); free(overallxmin);
        for (int i = 0; i < np - trans; i++) {
            free(poinc[i]);
        }
        free(poinc);
        free(xrms); free(overallxrms);
        free(customvalues);
    } // End of Parallel Block

    // Write results in file
    printf("\n\n  Writing Results in Output File...\n");
    EH_p_write_dyndiag_results(output_file, dim, nrms, rmsindex, results, pixels, ncustomvalues, customnames, nprintf, printfindex);
    // Free memory
    for (int i = 0; i < pixels; i++) {
        free(results[i]);
    }
    free(results);
    for (int i = 0; i < ncustomvalues; i++) {
            free(customnames[i]);
        }
    free(customnames);
}


void OS_full_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                        int indexX, int indexY, double *parrange, double *par, int npar, int nrms, int *rmsindex,
                                        void (*edosys)(int, double *, double, double *, double *),
                                        int ncustomvalues, int nprintf, int *printfindex,
                                        void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                                        int mode), int bifmode) {
    // Maximum length of custom names, if there is custom calculations
    size_t nchars = 20;
    // Declare matrix do store results
    int pixels = parrange[2]*parrange[5];  // Number of results
    double **results = malloc(pixels * sizeof **results);
    for (int i = 0; i < pixels; i++) {
        //results[i] = malloc((4 + (3*dim)) * sizeof **results);  // 4 for params, attrac, diffratac/ dim for xmax/ dim for xmin / dim for LE
        results[i] = malloc((4 + (5*dim) + (2*nrms) + nprintf) * sizeof **results);
    }
    // Declare rk4 timestep, final time, short initial time, pi and varstep
    double h, tf, s_T0;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control parameters
    double varstep[2];        
    varstep[0] = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    varstep[1] = (parrange[4] - parrange[3])/(parrange[5] - 1); // -1 in the denominator ensures the input resolution
    // Numerical control parameters
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations 
    // Declare variable to flag if all directions present same periodicity or not (0 = all the same, 1 = not the same)
    int diffAttrac = -1;
    // Declare variable to store attractor
    int attrac;
    // Prepare x vector to include perturbed values
    realloc_vector(x, ndim);
    // Allocate variable to store the names of the custom values
    char **customnames = malloc(ncustomvalues * sizeof *customnames);
    for (int i = 0; i < ncustomvalues; i++) {
        customnames[i] = malloc(nchars * sizeof *customnames);
    }
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, bifmode, ndiv, np, trans, ndim, maxper, varstep, npar, results, edosys, customfunc, customnames, nchars, pi, rmsindex, nrms, printfindex, nprintf, ncustomvalues) \
                                       private(t, h, tf, s_T0, attrac) \
                                       firstprivate(x, indexX, indexY, parrange, par, diffAttrac)
    {   
        //Get number of threads
        int ID = omp_get_thread_num();
        // Allocate memory for x` = f(x)
        double *f = malloc(dim * sizeof *f);
        // Allocate memory for vectors necessary for lyapunov exponents calculation
        double *cum = malloc(dim * sizeof *cum);            // Cumulative Vector
        double *lambda = malloc(dim *sizeof *lambda);       // Lyapunov Exponents vector
        double *s_cum = malloc(dim * sizeof *s_cum);        // Short Cumulative Vector
        double *s_lambda = malloc(dim * sizeof *s_lambda);  // Short Lyapunov Exponents Vector
        double *znorm = malloc(dim * sizeof *znorm);        // Norm of Vectors
        double *gsc = malloc((dim - 1) * sizeof *gsc);      // Inner Products Vector
        // Declare vector and allocate memory to store poincare map values: poinc[number of permanent regime forcing periods][dimension original system]
        double **poinc = malloc((np - trans) * sizeof **poinc);
        for (int i = 0; i < np - trans; i++) {
            poinc[i] = malloc(dim * sizeof **poinc);
        }
        // Declare vector to store the chosen Lyapunov Exponents to determine the attractor
        double *LE = malloc(dim * sizeof *LE);
        // Declare vector for temporary storage of periodicity values to check if all directions are equal
        int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
        // Store Initial Conditions
        double t0 = t;
        double *IC = malloc(dim * sizeof *IC);
        for (int i = 0; i < dim; i++) {
            IC[i] = (*x)[i];
        }
        // Declare memory to store min and max values
        double *xmax = malloc(dim * sizeof *xmax);
        double *xmin = malloc(dim * sizeof *xmin);
        double *overallxmin = malloc(dim * sizeof *overallxmin);
        double *overallxmax = malloc(dim * sizeof *overallxmax);
        // Allocate RMS variables
        double *xrms = malloc(dim * sizeof *xrms);
        double *overallxrms = malloc(dim * sizeof *overallxrms);
        // Allocate variable to store custom values
        double *customvalues = malloc(ncustomvalues * sizeof *customvalues);
        // Mumber of integration steps
        int N = np*ndiv;
        // Percentage of integration steps considered steady state regime 
        double steadystateperc = 1 - ((double)trans/(double)np);
        // Convert function arguments as local (private) variables
        double *X = convert_argument_to_private(*x, ndim);
        double *PAR = convert_argument_to_private(par, npar);
        // Check if there is any custom names to be inserted on the header of the output file
        if (ncustomvalues > 0) {
            customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 0);
        }
        // Index to identify position to write results
        int index;                                          
        // Declare counter for parallelized loop
        int k, m;
        #pragma omp for schedule(static) private(k, m)
        // Starts the parallel loop for Y control parameter
        for (k = 0; k < (int)parrange[5]; k++) {
            PAR[indexY] = parrange[3] + k*varstep[1]; // Increment value
            // Reset Initial conditions for the beggining of a horizontal line
            for (int i = 0; i < dim; i++) {
                X[i] = IC[i];
            }
            // Starts the loop for X control parameter
            for (m = 0; m < (int)parrange[2]; m++) {
                PAR[indexX] = parrange[0] + m*varstep[0]; // Increment Value
                // Check the mode of the diagram
                if (bifmode == 1) {
                    // Reset Initial conditions in each diagram step
                    for (int i = 0; i < dim; i++) {
                        X[i] = IC[i];
                    }
                }
                // Reset Variables
                t = t0;
                for (int i = 0; i < dim; i++) {
                    lambda[i] = 0.0;
                    s_lambda[i] = 0.0;
                    LE[i] = 0.0;
                    xrms[i] = 0.0;
                    overallxrms[i] = 0.0;
                    xmin[i] = 0.0;
                    overallxmin[i] = X[i];
                    xmax[i] = 0.0;
                    overallxmax[i] = X[i];
                }
                if (ncustomvalues > 0) {
                    for (int i = 0; i < ncustomvalues; i++) {
                        customvalues[i] = 0.0;
                    }
                }
                // Vary timestep if varpar = par[0], varying also final time and short initial time
                h = (2 * pi) / (ndiv * PAR[0]);              // par[0] = OMEGA
                tf = h*np*ndiv;                              // Final time
                s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time
                // Assign initial perturbation
                perturb_wolf(&X, dim, ndim, &cum, &s_cum);
                // Call Runge-Kutta 4th order integrator N = nP * nDiv times
                for (int i = 0; i < np; i++) {
                    for (int j = 0; j < ndiv; j++) {
                        rk4(ndim, X, t, h, PAR, f, edosys);
                        lyapunov_wolf(&X, t, h, dim, ndim, s_T0, &cum, &s_cum, &lambda, &s_lambda, &znorm, &gsc);
                        t = t + h;
                        // Apply poincare map at permanent regime
                        if (i >= trans) {
                            // Get max and min values at permanent regime
                            for (int q = 0; q < dim; q++) {
                                // Initialize xmax[dim] and xmin[dim] with first values of x[dim] at permanent regime
                                if (i == trans && j == 0) {
                                    xmax[q] = X[q];
                                    xmin[q] = X[q];
                                }
                                max_value(X[q], &xmax[q]);
                                min_value(X[q], &xmin[q]);
                            }
                            // Accumulate squared values to RMS computation in permanent regime
                            if (nrms > 0) {
                                for (int q = 0; q < nrms; q++) {
                                    xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], X[rmsindex[q]], N, 0);
                                }
                            }
                            // Perform custom calculations in steady state regime (if there is calculations to be done)
                            if (ncustomvalues > 0) {
                                customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 1);
                            }
                            // Choose any point in the trajectory for poincare section placement
                            if (j == 1) {
                                // Stores poincare values in poinc[np - trans][dim] vector
                                for (int p = 0; p < dim; p++) {
                                    poinc[i - trans][p] = X[p];
                                }
                            }
                        }
                        // Get overall max and min values
                        for (int q = 0; q < dim; q++) {
                            max_value(X[q], &overallxmax[q]);
                            min_value(X[q], &overallxmin[q]);
                        }
                        // Accumulate squared values for RMS computation for all time domain
                        if (nrms > 0) {
                            for (int q = 0; q < nrms; q++) {
                                overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], X[rmsindex[q]], N, 0);
                            }
                        }
                        // Perform custom calculations over the entire time series (transient + steady state), if there is calculations to be done
                        if (ncustomvalues > 0) {
                            customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 2);
                        }
                    }
                }
                // Compute RMS values of state variables
                if (nrms > 0) {
                    for (int q = 0; q < nrms; q++) {
                        xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], X[rmsindex[q]], N, 1);
                        overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], X[rmsindex[q]], N, 1);
                    }
                }
                // Perform custom calculations at the end of the time series (if there is calculations to be done)
                if (ncustomvalues > 0) {
                    customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 3);
                }
                // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
                store_LE(dim, lambda, s_lambda, LE);
                // Verify the type of motion of the system
                attrac = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
                // Write results in matrix
                index = (int)parrange[2]*k + m;
                results[index][0] = PAR[indexY];
                results[index][1] = PAR[indexX];
                results[index][2] = (double)attrac;
                results[index][3] = (double)diffAttrac;
                for (int r = 4; r < dim + 4; r++) {
                    results[index][r] = LE[r-4];
                }
                for (int r = dim + 4; r < 4 + (2*dim); r++) {
                    results[index][r] = xmax[r - 4 - dim];
                }
                for (int r = 4 + (2*dim); r < 4 + (3*dim); r++) {
                    results[index][r] = xmin[r - 4 - (2*dim)];
                }
                for (int r = 4 + (3*dim); r < 4 + (4*dim); r++) {
                    results[index][r] = overallxmax[r - 4 - (3*dim)];
                }
                for (int r = 4 + (4*dim); r < 4 + (5*dim); r++) {
                    results[index][r] = overallxmin[r - 4 - (4*dim)];
                }
                for (int r = 4 + (5*dim); r < 4 + (5*dim) + nrms; r++) {
                    results[index][r] = xrms[rmsindex[r - 4 - (5*dim)]];
                }
                for (int r = 4 + (5*dim) + nrms; r < 4 + (5*dim) + (2*nrms); r++) {
                    results[index][r] = overallxrms[rmsindex[r - 4 - (5*dim) - nrms]];
                }
                for (int r = 4 + (5*dim) + (2*nrms); r < 4 + (5*dim) + (2*nrms) + nprintf; r ++) {
                    results[index][r] = customvalues[printfindex[r - 4 - (5*dim) - (2*nrms)]];
                }
            }
            // Progress Monitor
            if (ID == 0) {
                //printf("PAR[%d] = %lf\n", indexX, PAR[indexX]);
                //printf("%-5s%-3d%-7s%-3d%-12s%-11lf%-11s%-12lf%-14s%-2d%-15s%-2d%-10s%-11lf%-10s%-11lf\n", "[k = ", k, "] [m = ", m, "]: parY = ", PAR[indexY], ", parX = ", PAR[indexX], ", Attractor = ", attrac, ", diffAttrac = ", diffAttrac, ", LE[0] = ", LE[0], ", LE[1] = ", LE[1]);
                progress_bar(0, PAR[indexY], parrange[3], (parrange[4] - varstep[1])/omp_get_num_threads());
                if (k == ((int)parrange[5] - 1)/omp_get_num_threads() ) {
                    progress_bar(1, PAR[indexY], parrange[3], (parrange[4] - varstep[1])/omp_get_num_threads());
                }
            }
        }
        // Free memory    
        free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
        free(znorm); free(gsc); free(LE); free(tmp_attrac); free(IC);
        free(xmax); free(xmin); free(overallxmax); free(overallxmin);
        for (int i = 0; i < np - trans; i++) {
            free(poinc[i]);
        }
        free(poinc);
        free(xrms); free(overallxrms);
        free(customvalues);
    } // End of Parallel Block
    
    // Write results in file
    printf("\n\n  Writing Results in Output File...\n");
    EH_p_write_fdyndiag_results(output_file, dim, nrms, rmsindex, results, pixels, ncustomvalues, customnames, nprintf, printfindex);
    // Free memory
    for (int i = 0; i < pixels; i++) {
        free(results[i]);
    }
    free(results);
    for (int i = 0; i < ncustomvalues; i++) {
            free(customnames[i]);
        }
    free(customnames);
}

void OS_full_forced_basin_of_attraction_2D_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                                    int indexX, int indexY, double *icrange, double *par, int npar, int nrms, int *rmsindex, 
                                                    void (*edosys)(int, double *, double, double *, double *), 
                                                    int ncustomvalues, int nprintf, int *printfindex,
                                                    void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                                                    int mode)) {
    // Maximum length of custom names, if there is custom calculations
    size_t nchars = 20;
    // Declare matrix do store results
    int pixels = icrange[2]*icrange[5];  // Number of results
    double **results = malloc(pixels * sizeof **results);
    for (int i = 0; i < pixels; i++) {
        results[i] = malloc((4 + (5*dim) + (2*nrms) + nprintf) * sizeof **results); 
    }
    // Declare rk4 timestep, final time, short initial time and pi 
    double h, tf, s_T0;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control initial conditions
    double icstep[2];        
    icstep[0] = (icrange[1] - icrange[0])/(icrange[2] - 1); // -1 in the denominator ensures the input resolution
    icstep[1] = (icrange[4] - icrange[3])/(icrange[5] - 1); // -1 in the denominator ensures the input resolution
    // Numerical control parameters
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations 
    // Declare variable to flag if all directions present same periodicity or not (0 = all the same, 1 = not the same)
    int diffAttrac = -1;
    // Declare variable to store attractor
    int attrac;
    // Prepare x vector to include perturbed values
    realloc_vector(x, ndim);
    // Allocate variable to store the names of the custom values
    char **customnames = malloc(ncustomvalues * sizeof *customnames);
    for (int i = 0; i < ncustomvalues; i++) {
        customnames[i] = malloc(nchars * sizeof *customnames);
    }
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, ndiv, np, trans, ndim, maxper, icstep, npar, results, edosys, customfunc, customnames, nchars, pi, rmsindex, nrms, printfindex, nprintf, ncustomvalues) \
                                       private(t, h, tf, s_T0, attrac) \
                                       firstprivate(x, indexX, indexY, icrange, par, diffAttrac)
    {   
        //Get number of threads
        int ID = omp_get_thread_num();
        // Allocate memory for x` = f(x)
        double *f = malloc(dim * sizeof *f);
        // Allocate memory for vectors necessary for lyapunov exponents calculation
        double *cum = malloc(dim * sizeof *cum);            // Cumulative Vector
        double *lambda = malloc(dim *sizeof *lambda);       // Lyapunov Exponents vector
        double *s_cum = malloc(dim * sizeof *s_cum);        // Short Cumulative Vector
        double *s_lambda = malloc(dim * sizeof *s_lambda);  // Short Lyapunov Exponents Vector
        double *znorm = malloc(dim * sizeof *znorm);        // Norm of Vectors
        double *gsc = malloc((dim - 1) * sizeof *gsc);      // Inner Products Vector
        // Declare vector and allocate memory to store poincare map values: poinc[number of permanent regime forcing periods][dimension original system]
        double **poinc = malloc((np - trans) * sizeof **poinc);
        for (int i = 0; i < np - trans; i++) {
            poinc[i] = malloc(dim * sizeof **poinc);
        }
        // Declare vector to store the chosen Lyapunov Exponents to determine the attractor
        double *LE = malloc(dim * sizeof *LE);
        // Declare vector for temporary storage of periodicity values to check if all directions are equal
        int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
        // Declare memory to store min and max values
        double *xmax = malloc(dim * sizeof *xmax);
        double *xmin = malloc(dim * sizeof *xmin);
        double *overallxmin = malloc(dim * sizeof *overallxmin);
        double *overallxmax = malloc(dim * sizeof *overallxmax);
        // Allocate RMS variables
        double *xrms = malloc(dim * sizeof *xrms);
        double *overallxrms = malloc(dim * sizeof *overallxrms);
        // Allocate variable to store custom values
        double *customvalues = malloc(ncustomvalues * sizeof *customvalues);
        // Mumber of integration steps
        int N = np*ndiv;
        // Percentage of integration steps considered steady state regime 
        double steadystateperc = 1 - ((double)trans/(double)np);
        // Allocate memory to store IC values
        double *IC = malloc(dim * sizeof *IC);
        // Convert function arguments as local (private) variables
        double *X = convert_argument_to_private(*x, ndim);
        double *PAR = convert_argument_to_private(par, npar);
        // Store Initial Conditions
        double t0 = t;
        for (int i = 0; i < dim; i++) {
            IC[i] = X[i];
        }
        // Check if there is any custom names to be inserted on the header of the output file
        if (ncustomvalues > 0) {
            customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 0);
        }
        // Index to identify position to write results
        int index;                                          
        // Declare counter for parallelized loop
        int k, m;
        #pragma omp for schedule(static) private(k, m)
        // Starts the parallel loop for Y control parameter
        for (k = 0; k < (int)icrange[5]; k++) {
            // Starts the loop for X control parameter
            for (m = 0; m < (int)icrange[2]; m++) {
                X[indexY] = icrange[3] + k*icstep[1]; // Increment value
                IC[indexY] = X[indexY];                    // Update IC value to write in result matrix
                X[indexX] = icrange[0] + m*icstep[0]; // Increment Value
                IC[indexX] = X[indexX];                    // Update IC value to write in result matrix
                // Reset Initial conditions in each basin step
                t = t0;
                for (int i = 0; i < dim; i++) {
                    X[i] = IC[i];
                }
                // Reset Variables
                for (int i = 0; i < dim; i++) {
                    lambda[i] = 0.0;
                    s_lambda[i] = 0.0;
                    LE[i] = 0.0;
                    xrms[i] = 0.0;
                    overallxrms[i] = 0.0;
                    xmin[i] = 0.0;
                    overallxmin[i] = X[i];
                    xmax[i] = 0.0;
                    overallxmax[i] = X[i];
                }
                if (ncustomvalues > 0) {
                    for (int i = 0; i < ncustomvalues; i++) {
                        customvalues[i] = 0.0;
                    }
                }
                // Vary timestep if varpar = par[0], varying also final time and short initial time
                h = (2 * pi) / (ndiv * PAR[0]);              // par[0] = OMEGA
                tf = h*np*ndiv;                              // Final time
                s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time
                // Assign initial perturbation
                perturb_wolf(&X, dim, ndim, &cum, &s_cum);
                // Call Runge-Kutta 4th order integrator N = nP * nDiv times
                for (int i = 0; i < np; i++) {
                    for (int j = 0; j < ndiv; j++) {
                        rk4(ndim, X, t, h, PAR, f, edosys);
                        lyapunov_wolf(&X, t, h, dim, ndim, s_T0, &cum, &s_cum, &lambda, &s_lambda, &znorm, &gsc);
                        t = t + h;
                        // Apply poincare map at permanent regime
                        if (i >= trans) {
                            // Get max and min values at permanent regime
                            for (int q = 0; q < dim; q++) {
                                // Initialize xmax[dim] and xmin[dim] with first values of x[dim] at permanent regime
                                if (i == trans && j == 0) {
                                    xmax[q] = X[q];
                                    xmin[q] = X[q];
                                }
                                max_value(X[q], &xmax[q]);
                                min_value(X[q], &xmin[q]);
                            }
                            // Accumulate squared values to RMS computation in permanent regime
                            if (nrms > 0) {
                                for (int q = 0; q < nrms; q++) {
                                    xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], X[rmsindex[q]], N, 0);
                                }
                            }
                            // Perform custom calculations in steady state regime (if there is calculations to be done)
                            if (ncustomvalues > 0) {
                                customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 1);
                            }
                            // Choose any point in the trajectory for poincare section placement
                            if (j == 1) {
                                // Stores poincare values in poinc[np - trans][dim] vector
                                for (int p = 0; p < dim; p++) {
                                    poinc[i - trans][p] = X[p];
                                }
                            }
                        }
                        // Get overall max and min values
                        for (int q = 0; q < dim; q++) {
                            max_value(X[q], &overallxmax[q]);
                            min_value(X[q], &overallxmin[q]);
                        }
                        // Accumulate squared values for RMS computation for all time domain
                        if (nrms > 0) {
                            for (int q = 0; q < nrms; q++) {
                                overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], X[rmsindex[q]], N, 0);
                            }
                        }
                        // Perform custom calculations over the entire time series (transient + steady state), if there is calculations to be done
                        if (ncustomvalues > 0) {
                            customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 2);
                        }
                    }
                }
                // Compute RMS values of state variables
                if (nrms > 0) {
                    for (int q = 0; q < nrms; q++) {
                        xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], X[rmsindex[q]], N, 1);
                        overallxrms[rmsindex[q]] = RMS(&overallxrms[rmsindex[q]], X[rmsindex[q]], N, 1);
                    }
                }
                // Perform custom calculations at the end of the time series (if there is calculations to be done)
                if (ncustomvalues > 0) {
                    customfunc(X, PAR, t, xrms, xmin, xmax, N, steadystateperc, ncustomvalues, customnames, nchars, customvalues, 3);
                }
                // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
                store_LE(dim, lambda, s_lambda, LE);
                // Verify the type of motion of the system
                attrac = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
                // Write results in matrix
                index = (int)icrange[2]*k + m;
                results[index][0] = IC[indexY];
                results[index][1] = IC[indexX];
                results[index][2] = (double)attrac;
                results[index][3] = (double)diffAttrac;
                for (int r = 4; r < dim + 4; r++) {
                    results[index][r] = LE[r-4];
                }
                for (int r = dim + 4; r < 4 + (2*dim); r++) {
                    results[index][r] = xmax[r - 4 - dim];
                }
                for (int r = 4 + (2*dim); r < 4 + (3*dim); r++) {
                    results[index][r] = xmin[r - 4 - (2*dim)];
                }
                for (int r = 4 + (3*dim); r < 4 + (4*dim); r++) {
                    results[index][r] = overallxmax[r - 4 - (3*dim)];
                }
                for (int r = 4 + (4*dim); r < 4 + (5*dim); r++) {
                    results[index][r] = overallxmin[r - 4 - (4*dim)];
                }
                for (int r = 4 + (5*dim); r < 4 + (5*dim) + nrms; r++) {
                    results[index][r] = xrms[rmsindex[r - 4 - (5*dim)]];
                }
                for (int r = 4 + (5*dim) + nrms; r < 4 + (5*dim) + (2*nrms); r++) {
                    results[index][r] = overallxrms[rmsindex[r - 4 - (5*dim) - nrms]];
                }
                for (int r = 4 + (5*dim) + (2*nrms); r < 4 + (5*dim) + (2*nrms) + nprintf; r ++) {
                    results[index][r] = customvalues[printfindex[r - 4 - (5*dim) - (2*nrms)]];
                }
                // Progress Monitor
                if (ID == 0) {
                    progress_bar(0, (double)k, 0, (icrange[5]-1)/omp_get_num_threads());
                    if (k == ((int)icrange[5] - 1)/omp_get_num_threads() ) {
                        progress_bar(1, (double)k, 0, (icrange[5]-1)/omp_get_num_threads());
                    }
                }
                //printf("results[%d][1] = IC[%d] | %lf = %lf\n", index, indexX, results[index][1], IC[indexX]);
            }
        }
        // Free memory    
        free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
        free(znorm); free(gsc); free(LE); free(tmp_attrac); free(IC);
        free(xmax); free(xmin); free(overallxmax); free(overallxmin);
        for (int i = 0; i < np - trans; i++) {
            free(poinc[i]);
        }
        free(poinc); 
        free(xrms); free(overallxrms);
        free(customvalues);
    } // End of Parallel Block
    
    // Write results in file
    printf("\n\n  Writing Results in Output File...\n");
    EH_p_write_fdyndiag_results(output_file, dim, nrms, rmsindex, results, pixels, ncustomvalues, customnames, nprintf, printfindex);

    // Free memory
    for (int i = 0; i < pixels; i++) {
        free(results[i]);
    }
    free(results);
    for (int i = 0; i < ncustomvalues; i++) {
            free(customnames[i]);
        }
    free(customnames);
}