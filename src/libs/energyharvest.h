#include <stdio.h>

double RMS(double *cum, double measure, int N, int mode);
void write_results_EH_ftimeseries(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, double Pout, int mode);
void EH_full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, double Pout, int mode));
