#include <stdio.h>

// Methods
double RMS(double *cum, double measure, int N, int mode);
void EH_print_RMS(FILE *info, int nRMS, int *rmsindex, double *xRMS, double *overallxRMS);

// Solutions
void EH_rk4_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, int nrms, int *rmsindex, double **xrms, double **overallxrms, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, int mode));
void EH_full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, int nrms, int *rmsindex, double **xrms, double **overallxrms, void (*edosys)(int, double *, double, double *, double *), void (*write_results)(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int mode));
void EH_bifurc_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, double t, double *x, int parindex, 
                     double *parrange, double *par, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                     void (*write_results)(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, int nrms, int *rmsindex, double *xrms, double *overallxrms, int mode), int bifmode);
void EH_full_bifurcation_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int maxper, double t,
                                  double **x, int parindex, double *parrange, double *par, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *),
                                  void (*write_results)(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *LE, int attractor, double **poinc, int diffattrac, int nrms, int *rmsindex, double *xrms, double *overallxrms, int mode), int bifmode);