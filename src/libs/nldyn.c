#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "edosystems.h"
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1 
#endif

// Methods
void rk4(int dim, double *x, double t, double h, double *par, double *f, void (*edosys)(int, double *, double, double *, double *)) {
    // Define arrays for temporary storage of the intermediary values of x[dim]
    double tmp[dim], k1[dim], k2[dim], k3[dim], k4[dim];
    // Calculate first slope
    edosys(dim, x, t, par, k1);
    // Assign next value for tmp[dim] to be inserted in the system of edos
    for (int i = 0; i < dim; i++) {
       tmp[i] = x[i] + 0.5 * (h * k1[i]);
    }
    // Calculate second slope
    edosys(dim, tmp, t + 0.5 * h, par, k2);
    // Assign next value for tmp[dim] to be inserted in the system of edos
    for (int i = 0; i < dim; i++) {
        tmp[i] = x[i] + 0.5 * (h * k2[i]);
    }
    // Calculate third slope
    edosys(dim, tmp, t + 0.5 * h, par, k3);
    // Assign next value for tmp[dim] to be inserted in the system of edos
    for (int i = 0; i < dim; i++) {
        tmp[i] = x[i] + h * k3[i];
    }
    // Calculate the fourth slope
    edosys(dim, tmp, t + h, par, k4);
    // Calculate the next value of x[dim]
    for (int i = 0; i < dim; i++) {
        x[i] = x[i] + (h/6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}

void realloc_vector(double **x, int ndim) {
    /* Realloc x[dim] vector to include new dimensions for linearized equations */
    double *tmp = realloc(*x, ndim * sizeof **x);
    if (tmp == NULL) {
        // Check if reallocation failed. If so, return..
        printf("Couldn't realloc x vector...\n");
        return;
    }
    else {
        (*x) = tmp;
    }
}

void perturb_wolf(double **x, int dim, int ndim, double **cum, double **s_cum) {
    /* Initial perturbation for linearized System of equations */
    // Attribute value of zero to x[dim] -> x[ndim] vetors
    for (int i = dim; i < ndim; i++) {
        (*x)[i] = 0.0;
    }
    for (int i = 0; i < dim; i ++) {
        // Perturb a single direction in the specific set of linearized equations
        (*x)[dim + ((i * dim) + i)] = 1.0;
        // Initialize cumulative ans short cumulative vectors with zero values
        (*cum)[i] = 0.0;
        (*s_cum)[i] = 0.0;
    }
    /* The user is responsible to free (x), (cum) and (s_cum) after the function call */
}

void lyapunov_wolf(double **x, double t, double h, int dim, int ndim, double s_t0, double **cum, double **s_cum, double **lambda, double **s_lambda, double **znorm, double **gsc) {
    // Define short time for short lyapunov
    double s_t = t - s_t0;
    /* Application of Gram-Schmidt Method to construct a new orthonormal basis */
    // Normalizes the 1st vector
    (*znorm)[0] = 0;
    for (int i = 0; i < dim; i++) {
        (*znorm)[0] = (*znorm)[0] + pow((*x)[dim + (i*dim)], 2);
    }
    (*znorm)[0] = sqrt((*znorm)[0]);
    
    for (int i = 0; i < dim; i++) {
        (*x)[dim + (i*dim)] = (*x)[dim + (i*dim)] / (*znorm)[0]; 
    }

    // Generates a new set of orthonormal vectors
    for (int i = 1; i < dim; i++) {
        // Generates (i-1) GSR coefficients (inner products)
        for (int j = 0; j < i; j++) {
            (*gsc)[j] = 0.0;

            for (int k = 0; k < dim; k++) {
                (*gsc)[j] = (*gsc)[j] + (*x)[dim + dim*k + i] * (*x)[dim + dim*k + j];
            }
        }
        // Build a new vector
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < i; k++) {
                (*x)[dim + dim*j + i] = (*x)[dim + dim*j + i] - (*gsc)[k]*(*x)[dim + dim*j + k];
            }
        }
        // Compute the norm of the new vector
        (*znorm)[i] = 0.0;

        for (int j = 0; j < dim; j++) {
            (*znorm)[i] = (*znorm)[i] + pow((*x)[dim + dim*j + i],2);
        }
        (*znorm)[i] = sqrt((*znorm)[i]);
        // Normalizes the new vector
        for (int j = 0; j < dim; j++) {
            (*x)[dim + dim*j + i] = (*x)[dim + dim*j + i] / (*znorm)[i];
        }
    }

    // Assigns magnitudes to current vectors, then calculates Lyapunov Exponents
    for (int i = 0; i < dim; i++) {
        (*cum)[i] = (*cum)[i] + (log((*znorm)[i]) / log(2));     // Cumulative vector 
        (*lambda)[i] = (*cum)[i] / (t + h);                   // Lyapunov Exponents
    }

    // Calculate Lyap Exponents in advance time 
    if (t >= s_t0) {
        for (int i = 0; i < dim; i++) {
            (*s_cum)[i] = (*s_cum)[i] + (log((*znorm)[i]) / log(2));     // Cumulative vector 
            (*s_lambda)[i] = (*s_cum)[i] / (s_t + h);                   // Lyapunov Exponents
        }
    }
}

int get_largest_element_in_array(int n, int *array) {
    // Safety check for pointer
    if (array == NULL) {
        printf("array is NULL, check input parameters or code...\n");
        exit(1);
    }
    // Create a temporary storage for the values of array
    int tmp_array[n];
    // Make tmp_array equals inserted array
    for (int i = 0; i < n; i++) {
        tmp_array[i] = array[i];
        //printf("tmp_array[%i] = %i\n", i, tmp_array[i]);
    }
    // Check the largest value in tmp_array
    for (int i = 1; i < n; i++) {
        if (tmp_array[0] < tmp_array[i]) {
            tmp_array[0] = tmp_array[i];
        }
    }
    // Returns the largest value of the array
    return tmp_array[0];
}

int check_periodicity(int dim, int np, double **poinc, int trans, int *tmp_attractor, int *diffattrac, int maxper) {
    // Add numerical tolerance
    double tol = 1e-4;                                          // Numerical tolerance
    // Check Parameters
    if (poinc == NULL || tmp_attractor == NULL) {
        printf("Poincare or/and Temporary Attractor vectors is/are NULL, check input parameters or code...\n");
        return -1;
    }
    // Check periodicity in each dimension
    for (int j = 0; j < dim; j++) {
        // Check periodicity of a specific j dimension
        for (int i = 1; i <= maxper; i++) {
            if ((fabs(poinc[np - trans - 1][j] - poinc[np - trans - 1 - i][j]) <= tol)) { 
                // If it finds a equal value, flag i
                tmp_attractor[j] = i;
                //printf("tmp_attractor[%i] = %i\n", j, i);
                // If the same occurs at 2*i, break and define periodicity
                if ((fabs(poinc[np - trans - 1][j] - poinc[np - trans - 1 - 2*i][j]) <= tol)) {
                    break;
                }
            } 
            // If i reaches the final number, attractor is period = maxper or greater
            if (i == maxper) {
                tmp_attractor[j] = i;
                //printf("tmp_attractor[%i] = %i\n", j, i);
            }
        }
    }
    // Check if all directions (dimensions) are periodic with the same periodicity
    for (int i = 0; i < dim - 1; i++) {
        for (int j = i + 1; j < dim; j++) {
            if (tmp_attractor[i] != tmp_attractor[j]) {
                // If not equal, flag it in the value of diffratac and get the largest periodicity
                tmp_attractor[0] = get_largest_element_in_array(dim, tmp_attractor);
                //printf("largest tmp_attractor element= %i\n", tmp_attractor[0]);
                (*diffattrac) = 1;
            }
            else {
                (*diffattrac) = 0;
            }
        } 
    }
    // return result of the first value of tmp_attractor (i.e.: the periodicity)
    return tmp_attractor[0];
}

int get_attractor(double **poinc, double *LE, int dim, int np, int trans, int *tmp_attractor, int *diffattrac, int maxper) {
    // Initialize variables
    int attractor;                                              // Returning value of the function
    int attracDiff;                                             // Variable that stores if periodicity is equal in all dimension (0 = equal, 1 = different)
    // Check Pointer
    if (LE == NULL) {
        printf("LE vectors is NULL, check input parameters or code...\n");
        return -1;
    }
    // Check if the system can show hyperchaotic motion or only chaotic
    if (dim <= 2) {
        // Check if system is periodic
        if (LE[0] < 0 ) {
            attractor = check_periodicity(dim, np, poinc, trans, tmp_attractor, diffattrac, maxper);
            //printf("attractor = %i\n", attractor);
        }
        // Check if the system is chaotic
        else if (LE[0] > 0) {
            attractor = maxper + 1;
        }
        // Error
        else {
            attractor = 0;          // Escape Point
            printf("Something went wrong trying to determine the type of motion of the system\n");
        }
    } 
    else if (dim > 2) {
        // Check if system is periodic
        if (LE[0] < 0 ) {
            attractor = check_periodicity(dim, np, poinc, trans, tmp_attractor, diffattrac, maxper);
            //printf("attractor = %i\n", attractor);
        }
        // Check if the system is chaotic
        else if (LE[0] > 0 && LE[1] <= 0) {
            attractor = maxper + 1;
        }
        // Check if the system is hyperchaotic
        else if (LE[0] > 0 && LE[1] > 0) {
            attractor = maxper + 2;
        }
        // Error
        else {
            attractor = 0;          // Escape Point
            printf("Something went wrong trying to determine the type of motion of the system\n");
        }
    }
    else {
        printf("Invalid dimension of the system...\n");
        return -1;
    }
    // Returns the value of the attractor
    return attractor;
}

void store_LE(int dim, double *lambda, double *s_lambda, double *le) {
    // Check if s_lambda is negative
    if (s_lambda[0] < 0) {
		for (int s = 0; s < dim; s++) {
			le[s] = s_lambda[s];
		}
	}
	else 
    {
		for (int s = 0; s < dim; s++) {
			le[s] = lambda[s];
		}
	}
}

double *convert_argument_to_private(double *arg, int nsize) {
    double *privarg = malloc(nsize * sizeof *privarg);
    for (int i = 0; i < nsize; i++) {
        privarg[i] = arg[i];
    }
    return privarg;
}

void add_one_row_2Dvector(double ***attrac, size_t *rows, size_t cols) {
    // Realloc matrix
    (*rows)++;
    double **tmp = realloc(*attrac, (*rows) * sizeof **attrac);
    if (tmp == NULL) {
        // Check if reallocation failed, If so, return...
        printf("Couldn't realloc 2D vector...\n");
        return;
    }
    else {
        for (int i = (*rows) - 1; i < (*rows); i++) {
            tmp[i] = malloc(cols * sizeof **attrac);
        }
        (*attrac) = tmp;
    }
}

bool check_if_array_is_all_one(int arr[], int dim) {
    for (int i = 0; i < dim; i++) {
        if (arr[i] != 1) {
            return false;
        }
    }
    return true;
}

void store_equilibrium_point(size_t *rows, size_t cols, double ***attrac, double *X, int dim, double tol) {
    // Define temporary variable to flag if attractor values are equal to X
    int flag[dim];
    // Check if there is no row in the attractor vector, if so, add a row with the first values of equilibrium positions
    if ((*rows) == 0) {
        add_one_row_2Dvector(attrac, rows, cols);
        for (int i = 0; i < dim; i++) {
            (*attrac)[0][i] = X[i];
        }
        (*attrac)[0][dim] = (*rows);
    } 
    else {
        // Check if values of X are equal of any value in attrac matrix. If so, return and don't add row to attrac
        for (int i = 0; i < (*rows); i++) {            
            for (int j = 0; j < dim; j ++) {
                if (fabs(X[j] - (*attrac)[i][j]) < tol) {
                    flag[j] = 1;
                } else {
                    flag[j] = 0;
                }
            }
            if (check_if_array_is_all_one(flag, dim) == true) {
                return;
            }         
        }
        // If there is no equal values, add row and assign the new equilibrium position
        add_one_row_2Dvector(attrac, rows, cols);
        for (int i = 0; i < dim; i++) {
            (*attrac)[(*rows) - 1][i] = X[i];
        }
        (*attrac)[(*rows) - 1][dim] = (*attrac)[(*rows) - 2][dim] + 1;   
    }
}

void max_value(double newvalue, double *oldvalue) {
    if (newvalue > (*oldvalue)) {
        (*oldvalue) = newvalue;
    }
}

void min_value(double newvalue, double *oldvalue) {
    if (newvalue < (*oldvalue)) {
        (*oldvalue) = newvalue;
    }
}

// Misc
void progress_bar(int mode, double var, double var_i, double var_f) {
    double perc;
    // Actual percentage
    if (mode == 1) {
        perc = 100;
    }
    else {
        perc = (var/(var_f - var_i))*100;
    }
    // Filled Part of the progress bar
    int fill = (perc * 50) / 100;  // 50 is the bar length
    printf("\r  Progress: |");
    for(int i = 0; i < fill; i++) {
        printf("#");
    }
    // Unfilled part of the progress bar
    for (int i = 0; i < 50 - fill; i++) {
        printf(" ");
    }
    if (perc > 100) {
        perc = 100;
    }
    printf("| %.2lf %%  (var = %.3lf)", perc, var);
    fflush(stdout);
}

void print_equilibrium_points(FILE* info, double **attrac, size_t rows, size_t cols, int dim) {
    printf("\n\n  Possible Stable Equilibrium Points: \n");
    fprintf(info, "\n  Possible Stable Equilibrium Points: \n");
    for (int j = 0; j < dim; j++) {
        if (j == 0) {
            printf("  ");
            fprintf(info, "  ");
        }
        printf("%-2s%-1d%-8s", "x[", j, "]");
        fprintf(info, "%-2s%-1d%-8s", "x[", j, "]");
    }
    printf("%-11s\n", "Attractor");
    fprintf(info, "%-11s\n", "Attractor");
    for (int i = 0; i < rows; ++i) { 
		for (int j = 0; j < cols; ++j) { 
            if (j == 0) {
                printf("  ");
                fprintf(info, "  ");
            }
            printf("%-10lf ", attrac[i][j]);
            fprintf(info, "%-10lf ", attrac[i][j]); 
        } 
		printf("\n");
        fprintf(info, "\n"); 
    }
}

// Solutions
void rk4_solution(FILE *output_file, int dim, int np, int ndiv, double t, double *x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, int mode)) {
    // Allocate x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Make the header o output file
    write_results(output_file, dim, t, x, 1);
    // Call Runge-Kutta 4th order integrator n = np * ndiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(dim, x, t, h, par, f, edosys);
            t = t + h;
            write_results(output_file, dim, t, x, 2);
        }
    }
    // Free Memory
    free(f);
}

void poincare_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, int mode)) {
    // Allocate x` = f(x)
    double *f = malloc(dim * sizeof *f);
    // Make the header of the output file
    write_results(output_file, dim, t, x, 0);
    // Call Runge-Kutta 4th order integrator n = np * ndiv times
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
                    write_results(output_file, dim, t, x, 2);
                }
            }
        }
    }
    // Free Memory
    free(f);
}

void lyap_wolf_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double **x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *lambda, double *s_lambda, int mode)) {
    // Allocate x` = f(x)
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
    double tf = h*np*ndiv;                              // Final time
    double s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time 
    // Prepare x vector to include perturbed values 
    realloc_vector(x, ndim);
    // Assign initial perturbation
    perturb_wolf(x, dim, ndim, &cum, &s_cum);
    /*printf("after perturb (outside function)\n");
    for(int i = 0; i < ndim; i++) {
        printf("x[%i] = %lf (memory: %p)\n", i, (*x)[i], &(*x)[i]);
    }*/
    // Make the header of output file and insert initial lyapunov exponents
    write_results(output_file, dim, t, lambda, s_lambda, 1);
    // Call Runge-Kutta 4th order integrator N = nP * nDiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(ndim, *x, t, h, par, f, edosys);
            lyapunov_wolf(x, t, h, dim, ndim, s_T0, &cum, &s_cum, &lambda, &s_lambda, &znorm, &gsc);
            t = t + h;
            write_results(output_file, dim, t, lambda, s_lambda, 2);
        }
    }
    // Free memory    
    free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
    free(znorm); free(gsc); 
}

void bifurc_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, double t, double *x, int parindex, 
                     double *parrange, double *par, void (*edosys)(int, double *, double, double *, double *), 
                     void (*write_results)(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, int mode), int bifmode) {
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
    // Declare timestep and pi
    double h;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control parameters
    double varstep;        
    varstep = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    // Make the header of the output file
    write_results(output_poinc_file, dim, par[parindex], x, xmin, xmax, 0);
    write_results(output_file, dim, par[parindex], x, xmin, xmax, 2);
    // Starts sweep the control parameter
    for (int k = 0; k < (int)parrange[2]; k++) {
        par[parindex] = parrange[0] + k*varstep; // Increment value
        // Reset Variables
        t = t0;
        // Check the mode of the bifurcation
        if (bifmode == 1) {
            // Reset Initial conditions in each bifurcation step
            for (int i = 0; i < dim; i++) {
                x[i] = IC[i];
            }
        }
        // Reset initial values to xmax and xmin based on initial values of x
        /*for (int i = 0; i < dim; i++) {
            xmax[i] = x[i];
            xmin[i] = x[i];
        }*/
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
                        write_results(output_poinc_file, dim, par[parindex], x, xmin, xmax, 1);
                    }
                }
            }
        }
        write_results(output_file, dim, par[parindex], x, xmin, xmax, 3);
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
    free(f); free(IC); free(xmax); free(xmin);
}

void full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int mode)) {
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
    /*printf("after perturb (outside function)\n");
    for(int i = 0; i < ndim; i++) {
        printf("x[%i] = %lf (memory: %p)\n", i, (*x)[i], &(*x)[i]);
    }*/
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
            }
        }
    }
    // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
    store_LE(dim, lambda, s_lambda, LE);
    for (int i = 0; i< dim; i++) {
        printf("  LE[%i] = %.10lf\n", i, LE[i]);
    }
    // Declare vector for temporary storage of periodicity values to check if all directions are equal
    int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
    // Declare variable to flag if all directions present same periodicity or not (0 = all the same, 1 = not the same)
    int diffAttrac = -1;
    // Verify the type of motion of the system
    (*attrac) = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
    printf("  diffAttrac = %i\n", diffAttrac);
    // Free memory    
    free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
    free(znorm); free(gsc); free(LE); free(tmp_attrac);
    for (int i = 0; i < np - trans; i++) {
        free(poinc[i]);
    }
    free(poinc);
}

void full_bifurcation_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int maxper, double t,
                               double **x, int parindex, double *parrange, double *par, void (*edosys)(int, double *, double, double *, double *),
                               void (*write_results)(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *LE, int attractor, double **poinc, int diffattrac, int mode), int bifmode) {
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
    write_results(output_poinc_file, dim, np, trans, par[parindex], (*x), xmin, xmax, LE, attrac, poinc, diffAttrac, 0);
    write_results(output_file, dim, np, trans, par[parindex], (*x), xmin, xmax, LE, attrac, poinc, diffAttrac, 2);
    // Starts to increment bifurcation control parameter
    for (int k = 0; k < (int)parrange[2]; k++) {
        par[parindex] = parrange[0] + k*varstep; // Increment value
        // Reset Variables
        t = t0;
        for (int i = 0; i < dim; i++) {
            lambda[i] = 0.0;
            s_lambda[i] = 0.0;
            LE[i] = 0.0;
        }
        // Check the mode of the bifurcation
        if (bifmode == 1) {
            // Reset Initial conditions in each bifurcation step
            for (int i = 0; i < dim; i++) {
                (*x)[i] = IC[i];
            }
        }
        // Reset initial values to xmax and xmin based on initial values of x
        /*for (int i = 0; i < dim; i++) {
            xmax[i] = (*x)[i];
            xmin[i] = (*x)[i];
        }*/
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
                    // Choose any point in the trajectory for poincare section placement
                    if (j == 1) {
                        // Stores poincare values in poinc[np - trans][dim] vector
                        for (int p = 0; p < dim; p++) {
                            poinc[i - trans][p] = (*x)[p];
                        }
                    }
                }
            }
        }
        // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
        store_LE(dim, lambda, s_lambda, LE);
        // Verify the type of motion of the system
        attrac = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
        // Write results in file
        write_results(output_poinc_file, dim, np, trans, par[parindex], (*x), xmin, xmax, LE, attrac, poinc, diffAttrac, 1);
        write_results(output_file, dim, np, trans, par[parindex], (*x), xmin, xmax, LE, attrac, poinc, diffAttrac, 3);
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
}

void dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x, int indexX, int indexY, double *parrange, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double varparX, double varparY, int attractor, double *LE, int diffattrac, int mode), int bifmode) {
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
    // Declare rk4 timestep, final time, short initial time, pi and varstep
    double h, tf, s_T0;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control parameters
    double varstep[2];        
    varstep[0] = (parrange[1] - parrange[0])/(parrange[2] - 1); // -1 in the denominator ensures the input resolution
    varstep[1] = (parrange[4] - parrange[3])/(parrange[5] - 1); // -1 in the denominator ensures the input resolution
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
    int attrac;
    // Declare vector for temporary storage of periodicity values to check if all directions are equal
    int *tmp_attrac = malloc(dim * sizeof *tmp_attrac);
    // Declare variable to flag if all directions present same periodicity or not (0 = all the same, 1 = not the same)
    int diffAttrac = -1;
    // Prepare x vector to include perturbed values
    realloc_vector(x, ndim);
    // Make the header of output files
    write_results(output_file, dim, par[indexX], par[indexY], attrac, LE, diffAttrac, 0);
    // Starts the loop for Y control parameter
    for (int k = 0; k < (int)parrange[5]; k++) {
        par[indexY] = parrange[3] + k*varstep[1]; // Increment value
        // Reset Initial conditions for the beggining of a horizontal line
        for (int i = 0; i < dim; i++) {
            (*x)[i] = IC[i];
        }
        // Starts the loop for X control parameter
        for (int m = 0; m < (int)parrange[2]; m++) {
            par[indexX] = parrange[0] + m*varstep[0]; // Increment Value
            // Reset Variables
            t = t0;
            for (int i = 0; i < dim; i++) {
                lambda[i] = 0.0;
                s_lambda[i] = 0.0;
                LE[i] = 0.0;
            }
            // Check the mode of the diagram
            if (bifmode == 1) {
                // Reset Initial conditions in each diagram step
                for (int i = 0; i < dim; i++) {
                    (*x)[i] = IC[i];
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
                        // Choose any point in the trajectory for poincare section placement
                        if (j == 1) {
                            // Stores poincare values in poinc[np - trans][dim] vector
                            for (int p = 0; p < dim; p++) {
                                poinc[i - trans][p] = (*x)[p];
                            }
                        }
                    }
                }
            }
            // Define which lyapunov will be taken: lambda[dim] or s_lambda[dim]
            store_LE(dim, lambda, s_lambda, LE);
            // Verify the type of motion of the system
            attrac = get_attractor(poinc, LE, dim, np, trans, tmp_attrac, &diffAttrac, maxper);
            // Write results in file
            write_results(output_file, dim, par[indexX], par[indexY], attrac, LE, diffAttrac, 1);
        }
        // Progress Monitor
        progress_bar(0, par[indexY], parrange[3], parrange[4]);
    }

    // Free memory    
    free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
    free(znorm); free(gsc); free(LE); free(tmp_attrac); free(IC);
    for (int i = 0; i < np - trans; i++) {
        free(poinc[i]);
    }
    free(poinc);
}

void parallel_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                         int indexX, int indexY, double *parrange, double *par, int npar,
                                         void (*edosys)(int, double *, double, double *, double *), int bifmode, 
                                         void (*write_results)(FILE *output_file, int dim, double **results, int pixels)) {
    // Declare matrix do store results
    int pixels = parrange[2]*parrange[5];  // Number of results
    double **results = malloc(pixels * sizeof **results);
    for (int i = 0; i < pixels; i++) {
        results[i] = malloc((4 + (3*dim)) * sizeof **results);  // 4 for params, attrac, diffratac/ dim for xmax/ dim for xmin / dim for LE
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
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, bifmode, ndiv, np, trans, ndim, maxper, varstep, npar, results, edosys, pi) \
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
        // Convert function arguments as local (private) variables
        double *X = convert_argument_to_private(*x, ndim);
        double *PAR = convert_argument_to_private(par, npar);
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
                // Reset Variables
                t = t0;
                for (int i = 0; i < dim; i++) {
                    lambda[i] = 0.0;
                    s_lambda[i] = 0.0;
                    LE[i] = 0.0;
                }
                // Check the mode of the diagram
                if (bifmode == 1) {
                    // Reset Initial conditions in each diagram step
                    for (int i = 0; i < dim; i++) {
                        X[i] = IC[i];
                    }
                }
                // Reset initial values to xmax and xmin based on initial values of x
                /*for (int i = 0; i < dim; i++) {
                    xmax[i] = X[i];
                    xmin[i] = X[i];
                }*/
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
                            // Choose any point in the trajectory for poincare section placement
                            if (j == 1) {
                                // Stores poincare values in poinc[np - trans][dim] vector
                                for (int p = 0; p < dim; p++) {
                                    poinc[i - trans][p] = X[p];
                                }
                            }
                        }
                    }
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
        free(xmax); free(xmin);
        for (int i = 0; i < np - trans; i++) {
            free(poinc[i]);
        }
        free(poinc);
    } // End of Parallel Block
    
    // Write results in file
    printf("\n\n  Writing Results in Output File...\n");
    write_results(output_file, dim, results, pixels);

    // Free memory
    for (int i = 0; i < pixels; i++) {
        free(results[i]);
    }
    free(results);
} 

void ep_basin_of_attraction_2D(FILE *output_file, FILE *info_file, int dim, int np, int ndiv, double t, double **x, int indexX, int indexY, double *icrange, double *par,
                               int npar, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, double **results, int pixels, int dim)) {
    // Declare matrix do store results
    int pixels = icrange[2]*icrange[5];  // Number of results
    double **results = malloc(pixels * sizeof **results);
    for (int i = 0; i < pixels; i++) {
        results[i] = malloc((3 + dim) * sizeof **results); 
    }
    // Declare rk4 timestep and pi
    double h;
    const double pi = 4 * atan(1);  // Pi number definition
    // Declare and define increment of control parameters
    double icstep[2];        
    icstep[0] = (icrange[1] - icrange[0])/(icrange[2] - 1); // -1 in the denominator ensures the input resolution
    icstep[1] = (icrange[4] - icrange[3])/(icrange[5] - 1); // -1 in the denominator ensures the input resolution
    // Declare variable to store attractors
    size_t rows = 0; size_t cols = dim + 1;
    double **attrac = malloc(rows * sizeof **attrac);
    for (int i = 0; i < rows; i++) {
        attrac[i] = malloc(cols * sizeof **attrac);
    }
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, ndiv, np, icstep, npar, results, edosys, attrac, rows, cols, pi) \
                                       private(t, h) \
                                       firstprivate(x, indexX, indexY, icrange, par)
    {   
        //Get number of threads
        int ID = omp_get_thread_num();
        // Allocate memory for x` = f(x)
        double *f = malloc(dim * sizeof *f);
        // Allocate memory to store IC values
        double *IC = malloc(dim * sizeof *IC);
        // Convert function arguments as local (private) variables
        double *X = convert_argument_to_private(*x, dim);
        double *PAR = convert_argument_to_private(par, npar);
        // Store Initial Conditions
        double t0 = t;
        for (int i = 0; i < dim; i++) {
            IC[i] = X[i];
        }
        // Index to identify position to write results
        int index;
        // Declare flag and tolerance to identify attractor
        int flag[dim]; 
        double tol = 1e-5;
        // Declare counter for parallelized loop
        int k, m;
        #pragma omp for schedule(static) private(k, m, flag) 
        // Starts the parallel loop for Y control parameter
        for (k = 0; k < (int)icrange[5]; k++) { 
            // Starts the loop for X control parameter
            for (m = 0; m < (int)icrange[2]; m++) {
                X[indexY] = icrange[3] + k*icstep[1];       // Increment value
                IC[indexY] = X[indexY];                     // Update IC value to write in result matrix
                X[indexX] = icrange[0] + m*icstep[0];       // Increment Value
                IC[indexX] = X[indexX];                     // Update IC value to write in result matrix                
                // Reset Variables
                t = t0;
                // Reset Initial conditions in each basin step
                for (int i = 0; i < dim; i++) {
                    X[i] = IC[i];
                }
                // Vary timestep if varpar = par[0], varying also final time and short initial time
                h = (2 * pi) / (ndiv * PAR[0]);              // par[0] = OMEGA
                // Call Runge-Kutta 4th order integrator N = nP * nDiv times
                for (int i = 0; i < np; i++) {
                    for (int j = 0; j < ndiv; j++) {
                        rk4(dim, X, t, h, PAR, f, edosys);
                        t = t + h;
                    }
                }
                // Compare X values to attrac vector and classify a new attractor if necessary
                // "critical" option to prevent multiple threads from accessing this section of code at the same time
                #pragma omp critical 
                {
                    store_equilibrium_point(&rows, cols, &attrac, X, dim, tol);
                }
                // Write corresponding results in matrix
                index = (int)icrange[2]*k + m;
                results[index][0] = IC[indexY];
                results[index][1] = IC[indexX];
                for (int i = 2; i < 2 + dim; i++) {
                    results[index][i] = X[i - 2];    
                }
                
                // Compare results of X to values in attractor matrix to classify the attractor
                for (int i = 0; i < rows; i++) {            
                    for (int j = 0; j < dim; j ++) {
                        if (fabs(X[j] - attrac[i][j]) < tol) {
                            flag[j] = 1;
                        } else {
                            flag[j] = 0;
                        }
                    }
                    if (check_if_array_is_all_one(flag, dim) == true) {
                        results[index][2 + dim] = attrac[i][dim];
                        break;
                    }         
                }
            }   
            // Progress Monitor
            if (ID == 0) {
                progress_bar(0, (double)k, 0, (icrange[5]-1)/omp_get_num_threads());
                if (k == ((int)icrange[5] - 1)/omp_get_num_threads() ) {
                    progress_bar(1, (double)k, 0, (icrange[5]-1)/omp_get_num_threads());
                }
            }
        }
        // Free memory    
        free(f); 
    } // End of Parallel Block
    
    //Print matrix containing the found attractors
    print_equilibrium_points(info_file, attrac, rows, cols, dim);
    // Write results in file
    printf("\n  Writing Results in Output File...\n");
    write_results(output_file, results, pixels, dim);

    // Free memory
    for (int i = 0; i < pixels; i++) {
        free(results[i]);
    }
    free(results);
}

void forced_basin_of_attraction_2D(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                         int indexX, int indexY, double *icrange, double *par, int npar,
                                         void (*edosys)(int, double *, double, double *, double *), 
                                         void (*write_results)(FILE *output_file, int dim, double **results, int pixels)) {
    // Declare matrix do store results
    int pixels = icrange[2]*icrange[5];  // Number of results
    double **results = malloc(pixels * sizeof **results);
    for (int i = 0; i < pixels; i++) {
        results[i] = malloc((4 + (3*dim)) * sizeof **results); 
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
    // Start of Parallel Block
    #pragma omp parallel default(none) shared(dim, ndiv, np, trans, ndim, maxper, icstep, npar, results, edosys, pi) \
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
                }
                // Reset initial values to xmax and xmin based on initial values of x
                /*for (int i = 0; i < dim; i++) {
                    xmax[i] = X[i];
                    xmin[i] = X[i];
                }*/
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
                            // Choose any point in the trajectory for poincare section placement
                            if (j == 1) {
                                // Stores poincare values in poinc[np - trans][dim] vector
                                for (int p = 0; p < dim; p++) {
                                    poinc[i - trans][p] = X[p];
                                }
                            }
                        }
                    }
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
                // Progress Monitor
                if (ID == 0) {
                    //printf("PAR[%d] = %lf\n", indexX, PAR[indexX]);
                    //printf("%-5s%-3d%-7s%-3d%-12s%-11lf%-11s%-12lf%-14s%-2d%-15s%-2d%-10s%-11lf%-10s%-11lf\n", "[k = ", k, "] [m = ", m, "]: parY = ", PAR[indexY], ", parX = ", PAR[indexX], ", Attractor = ", attrac, ", diffAttrac = ", diffAttrac, ", LE[0] = ", LE[0], ", LE[1] = ", LE[1]);
                    progress_bar(0, (double)k, 0, (icrange[5]-1)/omp_get_num_threads());
                    if (k == ((int)icrange[5] - 1)/omp_get_num_threads() ) {
                        progress_bar(1, (double)k, 0, (icrange[5]-1)/omp_get_num_threads());
                    }
                }
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
    } // End of Parallel Block
    
    // Write results in file
    printf("\n\n  Writing Results in Output File...\n");
    write_results(output_file, dim, results, pixels);

    // Free memory
    for (int i = 0; i < pixels; i++) {
        free(results[i]);
    }
    free(results);
}

// Not Finished (In Progress/tests)
void perturb_cldyn(double **x, int dim, int ndim, double perturb, double **cum, double **s_cum) {
    /* Initial perturbation for Cloned System of equations */
    // Attribute value of zero to x[dim] -> x[ndim] vetors
    for (int i = dim; i < ndim; i++) {
        (*x)[i] = 0.0;
    }
    for (int i = 0; i < dim; i ++) {
        // Perturb a single direction in the specific set of cloned equations
        (*x)[dim + ((i * dim) + i)] = (*x)[i] + perturb;
        // Initialize cumulative ans short cumulative vectors with zero values
        (*cum)[i] = 0.0;
        (*s_cum)[i] = 0.0;
    }
    /* The user is responsible to free (x), (cum) and (s_cum) after the function call */
}

void lyapunov_cldyn(double **x, double t, double h, int dim, int ndim, double perturb, double s_t0, double **cum, double **s_cum, double **lambda, double **s_lambda, double **znorm, double **gsc) {
    // Define short time for short lyapunov
    double s_t = t - s_t0;
    // Define the difference state vectors
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            (*x)[dim + dim*j + i] = (*x)[j] - (*x)[dim + dim*j + i];
            if(t < 2*h) {
                printf("x[%i] = x[%i] - x[%i] = %lf - %lf = %lf\n", dim + dim*j + i, j, dim + dim*j + i, (*x)[j], (*x)[dim + dim*j + i], (*x)[dim + dim*j + i]);
            }
        }
    }
    /* Application of Gram-Schmidt Method to construct a new orthonormal basis */
    // Normalizes the 1st vector
    (*znorm)[0] = 0;
    for (int i = 0; i < dim; i++) {
        (*znorm)[0] = (*znorm)[0] + pow((*x)[dim + (i*dim)], 2);
    }
    (*znorm)[0] = sqrt((*znorm)[0]);
    
    for (int i = 0; i < dim; i++) {
        (*x)[dim + (i*dim)] = (*x)[dim + (i*dim)] / (*znorm)[0]; 
    }

    // Generates a new set of orthonormal vectors
    for (int i = 1; i < dim; i++) {
        // Generates (i-1) GSR coefficients (inner products)
        for (int j = 0; j < i; j++) {
            (*gsc)[j] = 0.0;

            for (int k = 0; k < dim; k++) {
                (*gsc)[j] = (*gsc)[j] + (*x)[dim + dim*k + i] * (*x)[dim + dim*k + j];
            }
        }
        // Build a new vector
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < i; k++) {
                (*x)[dim + dim*j + i] = (*x)[dim + dim*j + i] - (*gsc)[k]*(*x)[dim + dim*j + k];
            }
        }
        // Compute the norm of the new vector
        (*znorm)[i] = 0.0;

        for (int j = 0; j < dim; j++) {
            (*znorm)[i] = (*znorm)[i] + pow((*x)[dim + dim*j + i],2);
        }
        (*znorm)[i] = sqrt((*znorm)[i]);
        // Normalizes the new vector
        for (int j = 0; j < dim; j++) {
            (*x)[dim + dim*j + i] = (*x)[dim + dim*j + i] / (*znorm)[i];
        }
    }

    // Assigns magnitudes to current vectors, then calculates Lyapunov Exponents in base 2 [bits/s]
    for (int i = 0; i < dim; i++) {
        (*cum)[i] = (*cum)[i] + (log((*znorm)[i]/perturb) / log(2));    // Cumulative vector 
        (*lambda)[i] = (*cum)[i] / (t + h);                             // Lyapunov Exponents
    }

    // Calculate Lyap Exponents in advance time in base 2 [bits/s]
    if (t >= s_t0) {
        for (int i = 0; i < dim; i++) {
            (*s_cum)[i] = (*s_cum)[i] + (log((*znorm)[i]/perturb) / log(2));    // Cumulative vector 
            (*s_lambda)[i] = (*s_cum)[i] / (s_t + h);                           // Lyapunov Exponents
        }
    }
    // Determine the new initial conditions in the orthogonal frame
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            (*x)[dim + dim*j + i] = (*x)[j] + (*x)[dim + dim*j + i]*perturb;
        }
    }
}

void lyap_cldyn_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double **x, double h, double *par, double perturb, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *lambda, double *s_lambda, int mode)) {
    // Allocate x` = f(x)
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
    double tf = h*np*ndiv;                              // Final time
    double s_T0 = ((double) trans/ (double) np) * tf;   // Advanced initial time 
    
    // Prepare x vector to include perturbed values and assign initial perturbation
    perturb_cldyn(x, dim, ndim, perturb, &cum, &s_cum);
    printf("after perturb (outside function)\n");
    for(int i = 0; i < ndim; i++) {
        printf("x[%i] = %lf (memory: %p)\n", i, (*x)[i], &(*x)[i]);
    }
    // Make the header of output file and insert initial lyapunov exponents
    write_results(output_file, dim, t, lambda, s_lambda, 1);
    // Call Runge-Kutta 4th order integrator N = nP * nDiv times
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ndiv; j++) {
            rk4(ndim, *x, t, h, par, f, edosys);
            lyapunov_cldyn(x, t, h, dim, ndim, perturb, s_T0, &cum, &s_cum, &lambda, &s_lambda, &znorm, &gsc);
            t = t + h;
            write_results(output_file, dim, t, lambda, s_lambda, 2);
        }
    }
    // Free memory    
    free(f); free(cum); free(s_cum); free(lambda); free(s_lambda);
    free(znorm); free(gsc); 
}


