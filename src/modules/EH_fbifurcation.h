#include <stdio.h>

void EH_fbifurc_read_params_and_IC(char *name, int *dim, int *npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *parindex, double **x, int *nrms, int **rmsindex, int *bifmode);
void EH_fbifurc_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int parindex, int nrms, int *rmsindex, int bifmode, char* edosys, char* mode);
void EH_fbifurcation(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));