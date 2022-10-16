#include <stdio.h>

void dyndiag(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));
void dyndiag_read_params_and_IC(char *name, int *dim, int *npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *indexX, int *indexY, double **x, int *bifmode);
void dyndiag_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int indexX, int indexY, int bifmode, char* funcname, char* mode);
