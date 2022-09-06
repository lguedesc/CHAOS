#include <stdio.h>

void EH_fts_read_params_and_IC(char *name, int *dim, int *npar, int *maxper, int *np, int *ndiv, int* trans, int *nrms, double *t, double **par, double **x, int **rmsindex);
void EH_fts_print_attractor(FILE* info, int attrac, int maxper, char *mode);
void EH_fts_print_info(FILE *info ,int dim, int npar, int maxper, int np, int ndiv, int trans, int nrms, double h, double t, double *x, double *par, int *rmsindex, char* funcname, char* mode);
void EH_print_RMS(FILE *info, int nRMS, int *rmsindex, double *xRMS, double *overallxRMS);
void EH_ftime_series(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));