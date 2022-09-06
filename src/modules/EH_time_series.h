#include <stdio.h>

void EH_tseries_read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, int *nrms, double *t, double **par, double **x, int **rmsindex);
void EH_tseries_print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, int nrms, double h, double t, double *x, double *par, int *rmsindex, char* funcname, char* mode);
void EH_timeseries(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));