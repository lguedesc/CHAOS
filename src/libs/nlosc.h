#ifndef NLOSC_H
#define NLOSC_H

#include <stdio.h>
#include "defines.h"

// Methods
double RMS(double *cum, double measure, int N, int mode);

// Functions to handle ang_info struct
ang_info *init_angle_struct(unsigned int nangles);
void free_ang_info_struct(ang_info *strct);

// Solutions
void HOS_timeseries_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, ang_info *angles, int nrms, int *rmsindex, double **xrms, double **overallxrms, 
                             double **xmin, double **xmax, double **overallxmin, double **overallxmax, void (*edosys)(int, double *, double, double *, double *),  
                             int ncustomvalues, char ***customnames, double **customvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex,
                             void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode));

void HOS_poincare_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, ang_info *angles, void (*edosys)(int, double *, double, double *, double *));

void HOS_full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, 
                                 ang_info *angles, int nrms, int *rmsindex, double **xrms, double **overallxrms, double **xmin, double **xmax, double **overallxmin, double **overallxmax,
                                 void (*edosys)(int, double *, double, double *, double *), 
                                 int ncustomvalues, char ***customnames, double **customvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex,
                                 void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode));

void HOS_bifurc_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, double t, double *x, int parindex, 
                        double *parrange, double *par, ang_info *angles, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                        int ncustomvalues, int nprintf, int *printfindex,
                        void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                        int mode), int bifmode);

void HOS_full_bifurcation_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int maxper, double t,
                                  double **x, int parindex, double *parrange, double *par, ang_info *angles, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                                  int ncustomvalues, int nprintf, int *printfindex,
                                  void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                                  int mode), int bifmode);

void HOS_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                    int indexX, int indexY, double *parrange, double *par, ang_info *angles, int npar, int nrms, int *rmsindex,
                                    void (*edosys)(int, double *, double, double *, double *),
                                    int ncustomvalues, int nprintf, int *printfindex,
                                    void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                                    int mode), int bifmode);

void HOS_full_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                        int indexX, int indexY, double *parrange, double *par, ang_info *angle, int npar, int nrms, int *rmsindex,
                                        void (*edosys)(int, double *, double, double *, double *),
                                        int ncustomvalues, int nprintf, int *printfindex,
                                        void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                                        int mode), int bifmode);

void HOS_full_forced_basin_of_attraction_2D_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                                    int indexX, int indexY, double *icrange, double *par, ang_info *angles, int npar, int nrms, int *rmsindex, 
                                                    void (*edosys)(int, double *, double, double *, double *), 
                                                    int ncustomvalues, int nprintf, int *printfindex,
                                                    void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue,
                                                    int mode));                        

#endif