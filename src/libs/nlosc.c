#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "iofiles.h"
#include "nldyn.h"
#include "odesolvers.h"
#include "basic.h"
#include "interface.h"
#include "defines.h"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1 
#endif

// General functions 
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

static void store_angle_results_in_matrix(double *vecmin, double *vecmax, ang_info *angles, int *col_offset, int index, double ***results) {
    // vecmax
    for (int r = 0; r < angles->n_angles; r++) {
        if (fabs(vecmax[angles->index[r]] - vecmin[angles->index[r]]) > TWOPI) {
            (*results)[index][(*col_offset) + r] = PI;
        }
        else {
            (*results)[index][(*col_offset) + r] = remainder(vecmax[angles->index[r]], TWOPI);
        }
    }
    (*col_offset) += angles->n_angles;
    // vecmin
    for (int r = 0; r < angles->n_angles; r++) {
        if (fabs(vecmax[angles->index[r]] - vecmin[angles->index[r]]) > TWOPI) {
            (*results)[index][(*col_offset) + r] = -PI;
        }
        else {
            (*results)[index][(*col_offset) + r] = remainder(vecmin[angles->index[r]], TWOPI);
        }
    }
    (*col_offset) += angles->n_angles;
}

static void store_results_in_matrix(double *vec, int vecsize, int *col_offset, int index, double ***results) {
    for (int r = 0; r < vecsize; r++) {
        (*results)[index][(*col_offset) + r] = vec[r];
    }
    (*col_offset) += vecsize;
}

static void store_specific_results_in_matrix(double *vec, int vecsize, int *vec_index, int *col_offset, int index, double ***results) {
    for (int r = 0; r < vecsize; r++) {
        (*results)[index][(*col_offset) + r] = vec[vec_index[r]];
    }
    (*col_offset) += vecsize;
}

static void store_results_in_matrix_dyndiag(double ***results, int k, int m, double *parrange, double *PAR, ang_info *angles, int indexY, int indexX,
                                            int attrac, int dim, double *xmax, double *xmin, double *overallxmax, double *overallxmin,
                                            int nrms, int *rmsindex, double *xrms, double *overallxrms, int nprintf, int *printfindex, 
                                            double *customvalues) {
    // Index to identify position in results array to store results
    int index = (int)parrange[2]*k + m;
    // Write Control parameters and attractor in matrix
    (*results)[index][0] = PAR[indexY];
    (*results)[index][1] = PAR[indexX];
    (*results)[index][2] = (double)attrac;
    int col_offset = 3;
    // Store xmax[dim] values
    store_results_in_matrix(xmax, dim, &col_offset, index, results);
    // Store xmin[dim] values
    store_results_in_matrix(xmin, dim, &col_offset, index, results);
    // Store overallxmax[dim] values
    store_results_in_matrix(overallxmax, dim, &col_offset, index, results);
    // Store overallxmin[dim] values
    store_results_in_matrix(overallxmin, dim, &col_offset, index, results);
    // Store xmax[angles->n_angles]_remainder and xmin[angles->n_angles]_remainder values of angles
    store_angle_results_in_matrix(xmin, xmax, angles, &col_offset, index, results);
    // Store overallxmax[angles->n_angles]_remainder and overallxmin[angles->n_angles]_remainder values of angles
    store_angle_results_in_matrix(overallxmin, overallxmax, angles, &col_offset, index, results);
    // Store xrms[nrms] values 
    store_specific_results_in_matrix(xrms, nrms, rmsindex, &col_offset, index, results);
    // Store overallxrms[nrms] values
    store_specific_results_in_matrix(overallxrms, nrms, rmsindex, &col_offset, index, results);
    // Store customvalues[nprintf] values
    store_specific_results_in_matrix(customvalues, nprintf, printfindex, &col_offset, index, results);
}

static void store_results_in_matrix_fdyndiag(double ***results, int k, int m, double *parrange, double *PAR, ang_info *angles, int indexY, int indexX,
                                             int attrac, int dim, double *LE, double *xmax, double *xmin, double *overallxmax, double *overallxmin,
                                             int nrms, int *rmsindex, double *xrms, double *overallxrms, int nprintf, int *printfindex, 
                                             double *customvalues) {
    // Index to identify position in results array to store results
    int index = (int)parrange[2]*k + m;
    // Write Control parameters and attractor in matrix
    (*results)[index][0] = PAR[indexY];
    (*results)[index][1] = PAR[indexX];
    (*results)[index][2] = (double)attrac;
    int col_offset = 3;
    // Store LE[dim] values
    store_results_in_matrix(LE, dim, &col_offset, index, results);
    // Store xmax[dim] values
    store_results_in_matrix(xmax, dim, &col_offset, index, results);
    // Store xmin[dim] values
    store_results_in_matrix(xmin, dim, &col_offset, index, results);
    // Store overallxmax[dim] values
    store_results_in_matrix(overallxmax, dim, &col_offset, index, results);
    // Store overallxmin[dim] values
    store_results_in_matrix(overallxmin, dim, &col_offset, index, results);
    // Store xmax[angles->n_angles]_remainder and xmin[angles->n_angles]_remainder values of angles
    store_angle_results_in_matrix(xmin, xmax, angles, &col_offset, index, results);
    // Store overallxmax[angles->n_angles]_remainder and overallxmin[angles->n_angles]_remainder values of angles
    store_angle_results_in_matrix(overallxmin, overallxmax, angles, &col_offset, index, results);
    // Store xrms[nrms] values 
    store_specific_results_in_matrix(xrms, nrms, rmsindex, &col_offset, index, results);
    // Store overallxrms[nrms] values
    store_specific_results_in_matrix(overallxrms, nrms, rmsindex, &col_offset, index, results);
    // Store customvalues[nprintf] values
    store_specific_results_in_matrix(customvalues, nprintf, printfindex, &col_offset, index, results);
}


// Functions to handle ang_info struct
ang_info *init_angle_struct(unsigned int nangles) {
    ang_info *info = malloc(sizeof(ang_info));
    info->n_angles = nangles;
    info->index = malloc(nangles * sizeof(*(info->index)));

    return info;
}

void free_ang_info_struct(ang_info *strct) {
    free(strct->index);
    free(strct);
}

// Solutions
void HOS_timeseries_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, ang_info *angles, int nrms, int *rmsindex, double **xrms, double **overallxrms, 
                             double **xmin, double **xmax, double **overallxmin, double **overallxmax, void (*edosys)(int, double *, double, double *, double *),  
                             int ncustomvalues, char ***customnames, double **customvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex,
                             void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode)) {
    // Allocate x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Allocate memory and store initial conditions
    double t0 = t;
    double *IC = malloc(dim * sizeof *IC);
    for (int i = 0; i < dim; i++) {
        IC[i] = x[i];
    }
    // Allocate RMS variables
    (*xrms) = malloc(dim * sizeof **xrms);
    (*overallxrms) = malloc(dim * sizeof **overallxrms);
    // Initialize xRMS[dim] and overallxRMS[dim] vector
    for (int i = 0; i < dim; i ++) {
        (*xrms)[i] = 0.0;
        (*overallxrms)[i] = 0.0;
    }
    // Allocate memory to store min, max, overall min and overall max values
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
        (*customnames) = alloc_string_array(ncustomvalues, MAX_CCALC_NAME_LEN);
    } 
    // Mumber of integration steps
    int N = np*ndiv;
    // Percentage of integration steps considered steady state regime 
    double steadystateperc = 1 - ((double)trans/(double)np);
    // Initialize currentstep
    int currenttimestep = 0;
    // Check if there is any custom names to be inserted on the header of the output file
    if (ncustomvalues > 0) {
        customfunc(x, par, t, (*xrms), (*xmin), (*xmax), IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, (*customnames), (*customvalues), 0);
    }
    // Make the header of the output file
    HOS_write_timeseries_results(output_file, dim, t, x, ncustomvalues, (*customnames), (*customvalues), nprintf, printfindex, angles, 1);
    // Call Runge-Kutta 4th order integrator n = np * ndiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            // Current timestep
            currenttimestep = ndiv*i +j;
            // Call integrator and add the timestep value to time variable
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
                    customfunc(x, par, t, (*xrms), (*xmin), (*xmax), IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, (*customnames), (*customvalues), 1);    
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
                customfunc(x, par, t, (*xrms), (*xmin), (*xmax), IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, (*customnames), (*customvalues), 2);
            }
            // Write results in output file
            HOS_write_timeseries_results(output_file, dim, t, x, ncustomvalues, (*customnames), (*customvalues), nprintf, printfindex, angles, 2);
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
        customfunc(x, par, t, (*xrms), (*xmin), (*xmax), IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, (*customnames), (*customvalues), 3);
    }
    // Free Memory
    free_mem(f, IC, NULL);
}

void HOS_poincare_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, ang_info *angles, void (*edosys)(int, double *, double, double *, double *)) {
    // Allocate x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Make the header of the output file
    HOS_write_poinc_results(output_file, dim, t, x, angles, 1);
    // Call Runge-Kutta 4th order integrator n = np * ndiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(dim, x, t, h, par, f, edosys);
            t = t + h;
            // Apply poincare map at permanent regime
            if (i > trans) {
                // Choose any point in the trajectory for plane placement
                if (j == 1) {
                    // Print the result in output file
                    HOS_write_poinc_results(output_file, dim, t, x, angles, 2);
                }
            }
        }
    }
    // Free Memory
    free(f);
}

void HOS_full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, 
                                 ang_info *angles, int nrms, int *rmsindex, double **xrms, double **overallxrms, double **xmin, double **xmax, double **overallxmin, double **overallxmax,
                                 void (*edosys)(int, double *, double, double *, double *), 
                                 int ncustomvalues, char ***customnames, double **customvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex,
                                 void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode)) {
    // Allocate memory for x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Allocate memory and store initial conditions
    double t0 = t;
    double *IC = malloc(dim * sizeof *IC);
    for (int i = 0; i < dim; i++) {
        IC[i] = (*x)[i];
    }
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
        (*customnames) = alloc_string_array(ncustomvalues, MAX_CCALC_NAME_LEN);
    }
    // Mumber of integration steps
    int N = np*ndiv;
    // Percentage of integration steps considered steady state regime 
    double steadystateperc = 1 - ((double)trans/(double)np);
    // Numerical control parameters
    double numtol = 1e-4;                               // Numerical tolerance for determine the attractor
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations
    double tf = h*np*ndiv;                              // Final time
    double s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time 
    // Declare vector and allocate memory to store poincare map values: poinc[number of permanent regime forcing periods][dimension original system]
    double **poinc = (double **)alloc_2D_array(np - trans, dim, sizeof(double));
     // Declare vector to store the chosen Lyapunov Exponents to determine the attractor
    double *LE = malloc(dim * sizeof *LE);
    // Prepare x vector to include perturbed values 
    realloc_vector(x, ndim);
    // Assign initial perturbation
    perturb_wolf(x, dim, ndim, &cum, &s_cum);
    // Initialize current time step variable
    int currenttimestep = 0;
    // Check if there is any custom names to be inserted on the header of the output file
    if (ncustomvalues > 0) {
        customfunc((*x), par, t, (*xrms), (*xmin), (*xmax), IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, (*customnames), (*customvalues), 0);
    }
    // Make the header of output files
    HOS_write_ftimeseries_results(output_ftimeseries_file, dim, t, (*x), lambda, s_lambda, ncustomvalues, (*customnames), (*customvalues), nprintf, printfindex, angles, 1);
    HOS_write_poinc_results(output_poinc_file, dim, t, (*x), angles, 1);
    // Call ODE solver N = nP * nDiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            currenttimestep = ndiv*i +j;
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
                    HOS_write_poinc_results(output_poinc_file, dim, t, (*x), angles, 2);
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
                    customfunc((*x), par, t, (*xrms), (*xmin), (*xmax), IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, (*customnames), (*customvalues), 1);    
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
                customfunc((*x), par, t, (*xrms), (*xmin), (*xmax), IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, (*customnames), (*customvalues), 2);
            }       
            // Write Results in output file
            HOS_write_ftimeseries_results(output_ftimeseries_file, dim, t, (*x), lambda, s_lambda, ncustomvalues, (*customnames), (*customvalues), nprintf, printfindex, angles, 2);
        }
    }
    // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
    store_LE(dim, lambda, s_lambda, LE);
    // Verify the type of motion of the system
    //(*attrac) = get_attractor_new(poinc, LE, dim, np, trans, maxper, (*xmin), (*xmax), numtol);
    (*attrac) = get_attractor(poinc, LE, dim, np, trans, maxper, (*xmin), (*xmax), numtol, angles);
    // Compute RMS values of state variables
    if (nrms > 0) {
        for (int q = 0; q < nrms; q++) {
            (*xrms)[rmsindex[q]] = RMS(&(*xrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
            (*overallxrms)[rmsindex[q]] = RMS(&(*overallxrms)[rmsindex[q]], (*x)[rmsindex[q]], N, 1);
        }
    }
    // Perform "end" type custom calculations if there is calculations to be done
    if (ncustomvalues > 0) {
        customfunc((*x), par, t, (*xrms), (*xmin), (*xmax), IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, (*customnames), (*customvalues), 3);
    }
    // Free memory 
    free_mem(f, cum, s_cum, lambda, s_lambda, znorm, gsc, LE, IC, NULL);
    free_2D_mem((void**)poinc, np - trans);   
}

void HOS_bifurc_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, double t, double *x, int parindex, 
                        double *parrange, double *par, ang_info *angles, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                        int ncustomvalues, int nprintf, int *printfindex,
                        void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                        int mode), int bifmode) {
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
    char **customnames = alloc_string_array(ncustomvalues, MAX_CCALC_NAME_LEN);
    // Mumber of integration steps
    int N = np*ndiv;
    // Percentage of integration steps considered steady state regime 
    double steadystateperc = 1 - ((double)trans/(double)np);
    // Declare timestep and pi
    double h;
    // Declare and define increment of control parameters
    double varstep;        
    varstep = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    // Initialize current time step variable
    int currenttimestep = 0;
    // Check if there is any custom names to be inserted on the header of the output file
    if (ncustomvalues > 0) {
        customfunc(x, par, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 0);
    }
    // Make the header of the output file
    HOS_write_bifurc_results(output_poinc_file, dim, par[parindex], x, xmin, xmax, overallxmin, overallxmax, nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, angles , 0);    
    HOS_write_bifurc_results(output_file, dim, par[parindex], x, xmin, xmax, overallxmin, overallxmax, nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, angles, 2);
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
        currenttimestep = 0;
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
        h = (2 * PI) / (ndiv * par[0]);         // par[0] = OMEGA
        // Call Runge-Kutta 4th order integrator n = np * ndiv times
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < ndiv; j++) {
                currenttimestep = ndiv*i +j;
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
                        HOS_write_bifurc_results(output_poinc_file, dim, par[parindex], x, xmin, xmax, overallxmin, overallxmax, nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, angles, 1);
                    }
                    // Accumulate squared values to RMS computation in permanent regime
                    if (nrms > 0) {
                        for (int q = 0; q < nrms; q++) {
                            xrms[rmsindex[q]] = RMS(&xrms[rmsindex[q]], x[rmsindex[q]], N, 0);
                        }
                    }
                    // Perform custom calculations in steady state regime (if there is calculations to be done)
                    if (ncustomvalues > 0) {
                        customfunc(x, par, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 1);    
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
                    customfunc(x, par, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 2);
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
            customfunc(x, par, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 3);
        }
        // Print results in output file
        HOS_write_bifurc_results(output_file, dim, par[parindex], x, xmin, xmax, overallxmin, overallxmax, nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, angles, 3);
        // Progress Monitor
        /*
        if (parrange[2] > 100) {
            if (k % 50 == 0) {
                progress_bar(0, par[parindex], parrange[0], parrange[1]);
            }
            if (k == parrange[2] - 1) {
                progress_bar(9, par[parindex], parrange[0], parrange[1]);
            }
        } else {
            progress_bar(0, par[parindex], parrange[0], parrange[1]);
        }*/
        if (k == parrange[2] - 1) {
            progress_bar(1, (double)k, 0.0, parrange[2]);
        }
        else {
            progress_bar(0, (double)k, 0.0, parrange[2]);
        }
    }
    // Free Memory
    free_mem(f, IC, xmax, xmin, xrms, overallxrms, overallxmin, overallxmax, customvalues, NULL);
    free_2D_mem((void**)customnames, ncustomvalues);
}

void HOS_full_bifurcation_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int maxper, double t,
                                  double **x, int parindex, double *parrange, double *par, ang_info *angles, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                                  int ncustomvalues, int nprintf, int *printfindex,
                                  void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                                  int mode), int bifmode) {
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
    char **customnames = alloc_string_array(ncustomvalues, MAX_CCALC_NAME_LEN);
    // Mumber of integration steps
    int N = np*ndiv;
    // Percentage of integration steps considered steady state regime 
    double steadystateperc = 1 - ((double)trans/(double)np);
    // Declare timestep, final time and short initial time
    double h, tf, s_T0;
    // Declare and define increment of control parameters
    double varstep;        
    varstep = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    // Numerical control parameters
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations 
    // Declare vector and allocate memory to store poincare map values: poinc[number of permanent regime forcing periods][dimension original system]
    size_t npoinc = np - trans;
    double **poinc = (double **)alloc_2D_array(npoinc, dim, sizeof(double));
    // Declare vector to store the chosen Lyapunov Exponents to determine the attractor
    double *LE = malloc(dim * sizeof *LE);
    // Declare variables related to the determination of the attractor
    int attrac = 0;
    double numtol = 1e-4;
    // Prepare x vector to include perturbed values
    realloc_vector(x, ndim);
    // Initialize current time step variable
    int currenttimestep = 0;
    // Check if there is any custom names to be inserted on the header of the output file
    if (ncustomvalues > 0) {
        customfunc((*x), par, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 0);
    }
    // Make the header of output files
    HOS_write_fbifurc_results(output_poinc_file, dim, np, trans, par[parindex], (*x), xmin, xmax, overallxmin, overallxmax, LE, attrac, npoinc, poinc,
                             nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, angles, 0);
    HOS_write_fbifurc_results(output_file, dim, np, trans, par[parindex], (*x), xmin, xmax, overallxmin, overallxmax, LE, attrac, npoinc, poinc,
                             nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, angles, 2);
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
        h = (2 * PI) / (ndiv * par[0]);              // par[0] = OMEGA
        tf = h*np*ndiv;                              // Final time
        s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time
        // Assign initial perturbation
        perturb_wolf(x, dim, ndim, &cum, &s_cum);
        // Call Runge-Kutta 4th order integrator N = nP * nDiv times
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < ndiv; j++) {
                currenttimestep = ndiv*i +j;
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
                        customfunc((*x), par, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 1);    
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
                    customfunc((*x), par, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 2);
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
            customfunc((*x), par, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 3);
        }
        // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
        store_LE(dim, lambda, s_lambda, LE);
        // Verify the type of motion of the system
        attrac = get_attractor(poinc, LE, dim, np, trans, maxper, xmin, xmax, numtol, angles);
        // Write results in file
        HOS_write_fbifurc_results(output_poinc_file, dim, np, trans, par[parindex], (*x), xmin, xmax, overallxmin, overallxmax, LE, attrac, npoinc, poinc,
                                  nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, angles, 1);
        HOS_write_fbifurc_results(output_file, dim, np, trans, par[parindex], (*x), xmin, xmax, overallxmin, overallxmax, LE, attrac, npoinc, poinc,
                                  nrms, rmsindex, xrms, overallxrms, ncustomvalues, customnames, customvalues, nprintf, printfindex, angles, 3);
        // Progress Monitor
        /*if (parrange[2] > 100) {
            if (k % 50 == 0) {
                progress_bar(0, par[parindex], parrange[0], parrange[1]);
            }
            if (k == parrange[2] - 1) {
                progress_bar(0, par[parindex], parrange[0], parrange[1]);
            }
        } else {
            progress_bar(0, par[parindex], parrange[0], parrange[1]);
        }*/
        if (k == parrange[2] - 1) {
            progress_bar(1, (double)k, 0.0, parrange[2]);
        }
        else {
            progress_bar(0, (double)k, 0.0, parrange[2]);
        }

    }
    // Free memory
    free_mem(f, cum, s_cum, lambda, s_lambda, znorm, gsc, LE, IC, xmax, xmin, overallxmax, overallxmin,
             xrms, overallxrms, NULL);    
    free_2D_mem((void**)poinc, npoinc);
    free_2D_mem((void**)customnames, ncustomvalues);
}

void HOS_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                    int indexX, int indexY, double *parrange, double *par, ang_info *angles, int npar, int nrms, int *rmsindex,
                                    void (*edosys)(int, double *, double, double *, double *),
                                    int ncustomvalues, int nprintf, int *printfindex,
                                    void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                                    int mode), int bifmode) {
    // Declare matrix do store results
    int pixels = parrange[2]*parrange[5];                                          // Number of results
    int ncols = (3 + (4*dim) + (4*angles->n_angles) + (2*nrms) + nprintf);         // Number of columns (3: CparX, CparY, attrac) + (4: xmin, xmax, overallxmin, overallxmax)*dim + (4: xmin_remainder, xmax_remainder, overallxmin_remainder, overallxmax_remainder)*dim + (2: xrms, overallxrms)*nrms + nprintf
    double **results = (double**)alloc_2D_array(pixels, ncols, sizeof(double));
    // Declare timestep and tolerance for attractor identification
    double h;
    double numtol = 1e-4;
    // Declare and define increment of control parameters
    double varstep[2];        
    varstep[0] = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    varstep[1] = (parrange[4] - parrange[3])/(parrange[5] - 1); // -1 in the denominator ensures the input resolution
    // Declare variable to store attractor
    int attrac;
    // Allocate variable to store the names of the custom values
    char **customnames = alloc_string_array(ncustomvalues, MAX_CCALC_NAME_LEN);
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, bifmode, ndiv, np, trans, maxper, varstep, npar, results, edosys, customfunc, customnames, rmsindex, nrms, printfindex, nprintf, ncustomvalues, numtol, angles) \
                                       private(h, attrac) \
                                       firstprivate(x, t, indexX, indexY, parrange, par)
    {   
        //Get number of threads
        int ID = omp_get_thread_num();
        int nThr = omp_get_num_threads();
        // Allocate memory for x` = f(x)
        double *f = malloc(dim * sizeof *f);
        // Declare vector and allocate memory to store poincare map values: poinc[number of permanent regime forcing periods][dimension original system]
        size_t npoinc = np - trans;
        double **poinc = (double **)alloc_2D_array(npoinc, dim, sizeof(double));
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
        // initialize current time step variable
        int currenttimestep = 0;
        // Check if there is any custom names to be inserted on the header of the output file
        if (ncustomvalues > 0) {
            customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 0);
        }
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
                h = (2 * PI) / (ndiv * PAR[0]);              // par[0] = OMEGA
                // Call Runge-Kutta 4th order integrator N = nP * nDiv times
                for (int i = 0; i < np; i++) {
                    for (int j = 0; j < ndiv; j++) {
                        currenttimestep = ndiv*i +j;
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
                                customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 1);
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
                            customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 2);
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
                    customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 3);
                }
                // Verify the type of motion of the system
                attrac = check_periodicity(dim, np, poinc, trans, maxper, xmin, xmax, numtol, angles);
                // Write results in matrix
                store_results_in_matrix_dyndiag(&results, k, m, parrange, PAR, angles, indexY, indexX, attrac, dim, xmax, xmin, overallxmax, overallxmin,
                                                nrms, rmsindex, xrms, overallxrms, nprintf, printfindex, customvalues);
            }
            // Progress Monitor
            if (ID == 0) {
                if (k == ((int)parrange[5] - 1)/nThr ) {
                    progress_bar(1, PAR[indexY], parrange[3], (parrange[4] - varstep[1])/nThr);
                }
                else {
                    progress_bar(0, PAR[indexY], parrange[3], (parrange[4] - varstep[1])/nThr);
                }
            }
        }
        // Free memory    
        free_mem(f, IC, xmin, xmax, overallxmin, overallxmax, xrms, overallxrms, customvalues, NULL);
        free_2D_mem((void **)poinc, npoinc);
    } // End of Parallel Block

    // Write results in file
    HOS_p_write_dyndiag_results(output_file, dim, nrms, rmsindex, results, pixels, ncustomvalues, customnames, nprintf, printfindex, angles);
    // Free memory
    free_2D_mem((void **)results, pixels);
    free_2D_mem((void **)customnames, ncustomvalues);
}

void HOS_full_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                        int indexX, int indexY, double *parrange, double *par, ang_info *angles, int npar, int nrms, int *rmsindex,
                                        void (*edosys)(int, double *, double, double *, double *),
                                        int ncustomvalues, int nprintf, int *printfindex,
                                        void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                                        int mode), int bifmode) {
    // Declare matrix do store results
    int pixels = parrange[2]*parrange[5];  // Number of results
    int ncols = (3 + (5*dim) + (4*angles->n_angles) + (2*nrms) + nprintf);         // Number of columns (3: CparX, CparY, attrac) + (5: xmin, xmax, overallxmin, overallxmax, LE)*dim + (2: xrms, overallxrms)*nrms + nprintf
    double **results = (double**)alloc_2D_array(pixels, ncols, sizeof(double));
    // Declare rk4 timestep, final time, short initial time and numtol
    double h, tf, s_T0;
    double numtol = 1e-4;
    // Declare and define increment of control parameters
    double varstep[2];        
    varstep[0] = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    varstep[1] = (parrange[4] - parrange[3])/(parrange[5] - 1); // -1 in the denominator ensures the input resolution
    // Numerical control parameters
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations 
    // Declare variable to store attractor
    int attrac;
    // Prepare x vector to include perturbed values
    realloc_vector(x, ndim);
    // Allocate variable to store the names of the custom values
    char **customnames = alloc_string_array(ncustomvalues, MAX_CCALC_NAME_LEN);
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, bifmode, ndiv, np, trans, ndim, maxper, varstep, npar, results, edosys, customfunc, customnames, rmsindex, nrms, printfindex, nprintf, ncustomvalues, numtol, angles) \
                                       private(h, tf, s_T0, attrac) \
                                       firstprivate(x, t, indexX, indexY, parrange, par)
    {   
        //Get number of threads
        int ID = omp_get_thread_num();
        int nThr = omp_get_num_threads();
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
        size_t npoinc = np - trans;
        double **poinc = (double **)alloc_2D_array(npoinc, dim, sizeof(double));
        // Declare vector to store the chosen Lyapunov Exponents to determine the attractor
        double *LE = malloc(dim * sizeof *LE);
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
        // Initialize current time step variable
        int currenttimestep = 0;
        // Check if there is any custom names to be inserted on the header of the output file
        if (ncustomvalues > 0) {
            customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 0);
        }
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
                h = (2 * PI) / (ndiv * PAR[0]);              // par[0] = OMEGA
                tf = h*np*ndiv;                              // Final time
                s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time
                // Assign initial perturbation
                perturb_wolf(&X, dim, ndim, &cum, &s_cum);
                // Call Runge-Kutta 4th order integrator N = nP * nDiv times
                for (int i = 0; i < np; i++) {
                    for (int j = 0; j < ndiv; j++) {
                        currenttimestep = ndiv*i + j;
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
                                customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 1);
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
                            customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 2);
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
                    customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 3);
                }
                // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
                store_LE(dim, lambda, s_lambda, LE);
                // Verify the type of motion of the system
                attrac = get_attractor(poinc, LE, dim, np, trans, maxper, xmin, xmax, numtol, angles);
                // Write results in matrix
                store_results_in_matrix_fdyndiag(&results, k, m, parrange, PAR, angles, indexY, indexX, attrac, dim, LE, xmax, xmin, overallxmax, overallxmin,
                                                 nrms, rmsindex, xrms, overallxrms, nprintf, printfindex, customvalues);
            }
            // Progress Monitor
            if (ID == 0) {
                //printf("PAR[%d] = %lf\n", indexX, PAR[indexX]);
                //printf("%-5s%-3d%-7s%-3d%-12s%-11lf%-11s%-12lf%-14s%-2d%-15s%-2d%-10s%-11lf%-10s%-11lf\n", "[k = ", k, "] [m = ", m, "]: parY = ", PAR[indexY], ", parX = ", PAR[indexX], ", Attractor = ", attrac, ", diffAttrac = ", diffAttrac, ", LE[0] = ", LE[0], ", LE[1] = ", LE[1]);
                if (k == ((int)parrange[5] - 1)/nThr) {
                    progress_bar(1, PAR[indexY], parrange[3], (parrange[4] - varstep[1])/nThr);
                } 
                else {
                    progress_bar(0, PAR[indexY], parrange[3], (parrange[4] - varstep[1])/nThr);
                }
            }
        }
        // Free memory    
        free_mem(f, cum, s_cum, lambda, s_lambda, znorm, gsc, LE, IC, xmax, xmin, overallxmax, overallxmin,
                  xrms, overallxrms, customvalues, NULL);
        free_2D_mem((void **)poinc, npoinc);
    } // End of Parallel Block
    
    // Write results in file
    HOS_p_write_fdyndiag_results(output_file, dim, nrms, rmsindex, results, pixels, ncustomvalues, customnames, nprintf, printfindex, angles);
    // Free memory
    free_2D_mem((void **)results, pixels);
    free_2D_mem((void **)customnames, ncustomvalues);
}

void HOS_full_forced_basin_of_attraction_2D_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                                    int indexX, int indexY, double *icrange, double *par, ang_info *angles, int npar, int nrms, int *rmsindex, 
                                                    void (*edosys)(int, double *, double, double *, double *), 
                                                    int ncustomvalues, int nprintf, int *printfindex,
                                                    void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                                                    int mode)) {
    // Declare matrix do store results
    int pixels = icrange[2]*icrange[5];              // Number of results
    int ncols = (3 + (5*dim) + (4*angles->n_angles) + (2*nrms) + nprintf);  // Number of columns (3: CparX, CparY, attrac) + (5: xmin, xmax, overallxmin, overallxmax, LE)*dim + (2: xrms, overallxrms)*nrms + nprintf
    double **results = (double**)alloc_2D_array(pixels, ncols, sizeof(double));
    // Declare rk4 timestep, final time, short initial time and numtol 
    double h, tf, s_T0;
    double numtol = 1e-4;
    // Declare and define increment of control initial conditions
    double icstep[2];        
    icstep[0] = (icrange[1] - icrange[0])/(icrange[2] - 1); // -1 in the denominator ensures the input resolution
    icstep[1] = (icrange[4] - icrange[3])/(icrange[5] - 1); // -1 in the denominator ensures the input resolution
    // Numerical control parameters
    int ndim = dim + (dim * dim);                       // Define new dimension to include linearized dynamical equations 
    // Declare variable to store attractor
    int attrac;
    // Prepare x vector to include perturbed values
    realloc_vector(x, ndim);
    // Allocate variable to store the names of the custom values
    char **customnames = alloc_string_array(ncustomvalues, MAX_CCALC_NAME_LEN);
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, ndiv, np, trans, ndim, maxper, icstep, npar, results, edosys, customfunc, customnames, rmsindex, nrms, printfindex, nprintf, ncustomvalues, numtol, angles) \
                                       private(h, tf, s_T0, attrac) \
                                       firstprivate(x, t, indexX, indexY, icrange, par)
    {   
        //Get number of threads
        int ID = omp_get_thread_num();
        int nThr = omp_get_num_threads();
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
        size_t npoinc = np - trans;
        double **poinc = (double **)alloc_2D_array(npoinc, dim, sizeof(double));
        // Declare vector to store the chosen Lyapunov Exponents to determine the attractor
        double *LE = malloc(dim * sizeof *LE);
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
        // Initialize current time step variable
        int currenttimestep = 0;
        // Check if there is any custom names to be inserted on the header of the output file
        if (ncustomvalues > 0) {
            customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 0);
        }
        // Declare counter for parallelized loop
        int k, m;
        #pragma omp for schedule(static) private(k, m)
        // Starts the parallel loop for Y control parameter
        for (k = 0; k < (int)icrange[5]; k++) {
            // Starts the loop for X control parameter
            for (m = 0; m < (int)icrange[2]; m++) {
                X[indexY] = icrange[3] + k*icstep[1];       // Increment value
                IC[indexY] = X[indexY];                     // Update IC value to write in result matrix
                X[indexX] = icrange[0] + m*icstep[0];       // Increment Value
                IC[indexX] = X[indexX];                     // Update IC value to write in result matrix
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
                h = (2 * PI) / (ndiv * PAR[0]);              // par[0] = OMEGA
                tf = h*np*ndiv;                              // Final time
                s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time
                // Assign initial perturbation
                perturb_wolf(&X, dim, ndim, &cum, &s_cum);
                // Call Runge-Kutta 4th order integrator N = nP * nDiv times
                for (int i = 0; i < np; i++) {
                    for (int j = 0; j < ndiv; j++) {
                        currenttimestep = ndiv*i +j;
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
                                customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 1);
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
                            customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 2);
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
                    customfunc(X, PAR, t, xrms, xmin, xmax, IC, t0, N, currenttimestep, steadystateperc, ncustomvalues, customnames, customvalues, 3);
                }
                // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
                store_LE(dim, lambda, s_lambda, LE);
                // Verify the type of motion of the system
                attrac = get_attractor(poinc, LE, dim, np, trans, maxper, xmin, xmax, numtol, angles);
                // Write results in matrix
                //store_results_in_matrix_fdyndiag(&results, k, m, icrange, IC, indexY, indexX, attrac, dim, LE, xmax, xmin, overallxmax, overallxmin,
                //                                 nrms, rmsindex, xrms, overallxrms, nprintf, printfindex, customvalues);
                // Progress Monitor
                if (ID == 0) {
                    if (k == ((int)icrange[5] - 1)/nThr) {
                        progress_bar(1, (double)k, 0, (icrange[5]-1)/nThr);
                    } 
                    else {
                        progress_bar(0, (double)k, 0, (icrange[5]-1)/nThr);
                    }
                }
                //printf("results[%d][1] = IC[%d] | %lf = %lf\n", index, indexX, results[index][1], IC[indexX]);
            }
        }
        // Free memory
        free_mem(f, cum, s_cum, lambda, s_lambda, znorm, gsc, LE, IC, xmax, xmin, overallxmax, overallxmin, xrms,
                 overallxrms, customvalues, NULL);    
        free_2D_mem((void **)poinc, npoinc);        
    } // End of Parallel Block
    
    // Write results in file
    HOS_p_write_fdyndiag_results(output_file, dim, nrms, rmsindex, results, pixels, ncustomvalues, customnames, nprintf, printfindex, angles);

    // Free memory
    free_2D_mem((void **)results, pixels);
    free_2D_mem((void **)customnames, ncustomvalues);
}
