#include <stdio.h>

void lyap_wolf_read_params_and_IC(char *name, int *dim, int *npar, int *np, int *ndiv, int *trans, double *t, double **par, double **x);
void lyap_wolf_print_info(FILE *info ,int dim, int npar, int np, int ndiv, int trans, double h, double t, double *x, double *par, char* funcname, char* mode);
void lyapunov_exp_wolf(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *));
