#ifndef NLDYN_H
#define NLDYN_H

#include <stdio.h>
#include <stdbool.h>

//Methods
void realloc_vector(double **x, int ndim);
void perturb_wolf(double **x, int dim, int ndim, double **cum, double **s_cum);
void lyapunov_wolf(double **x, double t, double h, int dim, int ndim, double s_t0, double **cum, double **s_cum, double **lambda, double **s_lambda, double **znorm, double **gsc);
void store_LE(int dim, double *lambda, double *s_lambda, double *le);
int get_largest_element_in_array(int n, int *array);
int check_periodicity(int dim, int np, double **poinc, int trans, int *tmp_attractor, int *diffattrac, int maxper);
int get_attractor(double **poinc, double *LE, int dim, int np, int trans, int *tmp_attractor, int *diffattrac, int maxper);
double *convert_argument_to_private(double *arg, int nsize);
void add_one_row_2Dvector(double ***attrac, size_t *rows, size_t cols);
bool check_if_array_is_all_one(int arr[], int dim);
void store_equilibrium_point(size_t *rows, size_t cols, double ***attrac, double *X, int dim, double tol);
void max_value(double newvalue, double *oldvalue);
void min_value(double newvalue, double *oldvalue);

double *get_system_tol(int dim, double *xmin, double *xmax);
int check_periodicity_new(int dim, int np, double **poinc, int trans, int maxper, double *xmin, double *xmax, double numtol);
int get_attractor_new(double **poinc, double *LE, int dim, int np, int trans, int maxper, double *xmin, double *xmax, double numtol);

// Misc
void print_equilibrium_points(FILE* info, double **attrac, size_t rows, size_t cols, int dim);
// Solutions
void timeseries_solution(FILE *output_file, int dim, int np, int ndiv, double t, double *x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, int mode));
void poincare_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, int mode));
void lyap_wolf_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double **x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *lambda, double *s_lambda, int mode));
void ep_basin_of_attraction_2D(FILE *output_file, FILE *info_file, int dim, int np, int ndiv, double t, double **x, int indexX, int indexY, double *icrange, double *par,
                               int npar, void (*edosys)(int, double *, double, double *, double *));
// In progress Functions (in development)
void convergence_test_solution(int dim, int N, double t, double tf, double *x, int ntries, double *par, 
                                void (*edosys)(int, double *, double, double *, double *));



// Not in use
/*
void full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int mode));
void bifurc_solution(FILE* output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, double t, double *x, int parindex, 
                     double *parrange, double *par, void (*edosys)(int, double *, double, double *, double *), 
                     void (*write_results)(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, int mode), int bifmode);

void full_bifurcation_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int maxper, double t,
                               double **x, int parindex, double *parrange, double *par, void (*edosys)(int, double *, double, double *, double *),
                               void (*write_results)(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *LE, int attractor, double **poinc, int diffattrac, int mode), int bifmode);

void dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x, int indexX, int indexY, double *parrange, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double varparX, double varparY, int attractor, double *LE, int diffattrac, int mode), int bifmode);
void parallel_full_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x, int indexX, int indexY, double *parrange, double *par, int npar, void (*edosys)(int, double *, double, double *, double *), int bifmode, void (*write_results)(FILE *output_file, int dim, double **results, int pixels));
void forced_basin_of_attraction_2D(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                         int indexX, int indexY, double *icrange, double *par, int npar,
                                         void (*edosys)(int, double *, double, double *, double *), 
                                         void (*write_results)(FILE *output_file, int dim, double **results, int pixels));

*/

// In Progress Functions (also not in use)
/*
void perturb_cldyn(double **x, int dim, int ndim, double perturb, double **cum, double **s_cum);
void lyapunov_cldyn(double **x, double t, double h, int dim, int ndim, double perturb, double s_t0, double **cum, double **s_cum, double **lambda, double **s_lambda, double **znorm, double **gsc);
void lyap_cldyn_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double **x, double h, double *par, double perturb, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *lambda, double *s_lambda, int mode));

void parallel_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                         int indexX, int indexY, double *parrange, double *par, int npar,
                                         void (*edosys)(int, double *, double, double *, double *), int bifmode, 
                                         void (*write_results)(FILE *output_file, int dim, double **results, int pixels));
*/


#endif