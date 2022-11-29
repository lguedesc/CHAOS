#include <stdio.h>

// Methods
double RMS(double *cum, double measure, int N, int mode);

// Solutions
void OS_rk4_solution(FILE *output_file, int dim, int np, int ndiv, int trans, double t, double *x, double h, double *par, int nrms, int *rmsindex, double **xrms, double **overallxrms, 
                     double **xmin, double **xmax, double **overallxmin, double **overallxmax, void (*edosys)(int, double *, double, double *, double *),  
                     int ncustomvalues, char ***customnames, double **customvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex,
                     void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode));

void OS_full_timeseries_solution(FILE *output_ftimeseries_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int *attrac, int maxper, double t, double **x, double h, double *par, 
                                 int nrms, int *rmsindex, double **xrms, double **overallxrms, double **xmin, double **xmax, double **overallxmin, double **overallxmax,
                                 void (*edosys)(int, double *, double, double *, double *), 
                                 int ncustomvalues, char ***customnames, double **customvalues, int nprintf, int *printfindex, int nprintscr, int *printscrindex,
                                 void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode));

void OS_bifurc_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, double t, double *x, int parindex, 
                        double *parrange, double *par, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                        int ncustomvalues, int nprintf, int *printfindex,
                        void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                        int mode), int bifmode);

void OS_full_bifurcation_solution(FILE *output_file, FILE *output_poinc_file, int dim, int np, int ndiv, int trans, int maxper, double t,
                                  double **x, int parindex, double *parrange, double *par, int nrms, int *rmsindex, void (*edosys)(int, double *, double, double *, double *), 
                                  int ncustomvalues, int nprintf, int *printfindex,
                                  void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                                  int mode), int bifmode);

void OS_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                    int indexX, int indexY, double *parrange, double *par, int npar, int nrms, int *rmsindex,
                                    void (*edosys)(int, double *, double, double *, double *),
                                    int ncustomvalues, int nprintf, int *printfindex,
                                    void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                                    int mode), int bifmode);

void OS_full_dynamical_diagram_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                        int indexX, int indexY, double *parrange, double *par, int npar, int nrms, int *rmsindex,
                                        void (*edosys)(int, double *, double, double *, double *),
                                        int ncustomvalues, int nprintf, int *printfindex,
                                        void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                                        int mode), int bifmode);

void OS_full_forced_basin_of_attraction_2D_solution(FILE *output_file, int dim, int np, int ndiv, int trans, int maxper, double t, double **x,
                                                    int indexX, int indexY, double *icrange, double *par, int npar, int nrms, int *rmsindex, 
                                                    void (*edosys)(int, double *, double, double *, double *), 
                                                    int ncustomvalues, int nprintf, int *printfindex,
                                                    void (*customfunc)(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue,
                                                    int mode));                        
