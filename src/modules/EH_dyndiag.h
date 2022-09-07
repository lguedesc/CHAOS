#include <stdio.h>

void EH_dyndiag_read_params_and_IC(char *name, int *dim, int *npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *indexX, int *indexY, double **x, int *nrms, int **rmsindex, int *bifmode);
void EH_dyndiag_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int indexX, int indexY, int nrms, int *rmsindex, int bifmode, char* funcname, char* mode);
void EH_dyndiag(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));