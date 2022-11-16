#include <stdio.h>
#include <stdbool.h>

// Handles the reading, creation and checking of files and directories
char *convert_dir(const char *rawdir);
void OSmkdir(char *path);
void create_dir(const char *pathname);
void directory_exists(const char *pathname);
bool file_exists(const char* filename);
FILE *create_output_file(char *name, const char *ext, const char *dir);
char *get_input_filename(void);

// Write Results in file
void write_results(FILE *output_file, int dim, double t, double *x, int mode);
void write_results_lyap(FILE *output_file, int dim, double t, double *lambda, double *s_lambda, int mode);
void write_results_ftimeseries(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int mode);
//void write_bifurc_results(FILE *output_file, int dim, double varpar, double *x, int mode);
void write_bifurc_results(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, int mode);
//void write_fbifurc_results(FILE *output_file, int dim, int np, int trans, double varpar, double *x, int attractor, double **poinc, int diffattrac, int mode);
void write_fbifurc_results(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *LE, int attractor, double **poinc, int diffattrac, int mode);
void write_dyndiag_results(FILE *output_file, int dim, double varparX, double varparY, int attractor, double *LE, int diffattrac, int mode);
void p_write_fdyndiag_results(FILE *output_file, int dim, double **results, int pixels);
void p_write_dyndiag_results(FILE *output_file, int dim, double **results, int pixels);
void p_write_epbasin_results(FILE *output_file, double **results, int pixels, int dim);

// Write results in file for Energy Harvesting Toolbox
void EH_write_bifurc_results(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int nrms, int *rmsindex, double *xrms, double *overallxrms,
                             int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, int mode);
void EH_write_fbifurc_results(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *overallxmin, double *overallxmax, double *LE, int attractor, double **poinc, int diffattrac, int nrms, int *rmsindex, double *xrms, double *overallxrms,
                              int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, int mode);
void EH_p_write_fdyndiag_results(FILE *output_file, int dim, int nrms, int *rmsindex, double **results, int pixels);
void EH_p_write_dyndiag_results(FILE *output_file, int dim, int nrms, int *rmsindex, double **results, int pixels, int ncustomvalues, char **customnames, int nprintf, int *printfindex);
void EH_write_timeseries_results(FILE *output_file, int dim, double t, double *x, int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, int mode);
void EH_write_ftimeseries_results(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, int mode);