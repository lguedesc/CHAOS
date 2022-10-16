#include <stdio.h>

void fdyndiag_read_params_and_IC(char *name, int *dim, int *npar, int *maxper,  int *np, int *ndiv, int *trans, double *t, double **par, double **parrange, int *indexX, int *indexY, double **x, int *bifmode);
void fdyndiag_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, double t, double *x, double *par, double *parrange, int indexX, int indexY, int bifmode, char* edosys, char* mode);
void fdyndiag(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));