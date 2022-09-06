#include <stdio.h>

void tseries_read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, double *t, double **par, double **x);
void tseries_print_info(FILE *info ,int dim, int npar, int np, int ndiv, double h, double t, double *x, double *par, char *edosys, char *mode);
void timeseries(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));