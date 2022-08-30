#include <stdio.h>

void epbasin_read_params_and_IC(char *name, int *dim, int *npar,  int *np, int *ndiv, double *t, double **par, double **icrange, int *indexX, int *indexY, double **x);
void epbasin_print_info(FILE *info ,int dim, int npar, int np, int ndiv, double t, double *x, double *par, double *icrange, int indexX, int indexY, char* funcname, char* mode);
void epbasin(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));
