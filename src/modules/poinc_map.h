#include <stdio.h>

void poinc_read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, double *t, double **par, double **x);
void poinc_print_info(FILE *info ,int dim, int npar, int np, int ndiv, double h, double t, double *x, double *par, char *edosys, char *mode);
void poincaremap(char *funcname, void (*edosys)(int, double *, double, double *, double *));