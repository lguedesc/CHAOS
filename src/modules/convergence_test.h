#include <stdio.h>

void convergence_test(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));
void convergence_read_params_and_IC(char *name, int *dim, int *npar, int *n, double *t, double *tf, int *ntries, double **par, double **x);
void convergence_print_info(FILE *info ,int dim, int npar, int np, int ndiv, double h, double t, double *x, double *par, char* funcname, char* mode);