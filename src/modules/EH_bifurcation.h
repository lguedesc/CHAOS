#include <stdio.h>

void EH_bifurc_read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *parindex, double **x, int *nrms, int **rmsindex, int *bifmode);
void EH_bifurc_print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int parindex, int nrms, int *rmsindex, int bifmode, char* funcname, char* mode);
void EH_bifurcation(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));