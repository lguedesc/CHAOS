#include <stdio.h>

void fbifurc_read_params_and_IC(char *name, int *dim, int *npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *parindex, double **x, int *bifmode);
void fbifurc_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, double t, double *x, double *par, double *parrange, int parindex, int bifmode, char* edosys, char* mode);
void fbifurcation(char *funcname, void (*edosys)(int, double *, double, double *, double *));
