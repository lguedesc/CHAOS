#include <stdio.h>

void bifurc_read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *parindex, double **x, int *bifmode);
void bifurc_print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int parindex, int bifmode, char* funcname, char* mode);
void bifurcation(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));

