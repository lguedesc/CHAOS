#include <stdio.h>

void EH_forcedbasin_read_params_and_IC(char *name, int *dim, int *npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **icrange, int *indexX, int *indexY, double **x, int *nrms, int **rmsindex);
void EH_forcedbasin_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, double t, double *x, double *par, double *icrange, int indexX, int indexY, int nrms, int *rmsindex, char* funcname, char* mode);
void EH_forcedbasin(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));
