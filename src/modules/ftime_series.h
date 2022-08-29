#include <stdio.h>

void fts_read_params_and_IC(char *name, int *dim, int *npar, int *maxper, int *np, int *ndiv, int* trans, double *t, double **par, double **x);
void fts_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, double h, double t, double *x, double *par, char* edosys, char* mode);
void fts_print_attractor(FILE* info, int attrac, int maxper, char *mode);
void ftime_series(char *funcname, void (*edosys)(int, double *, double, double *, double *));
