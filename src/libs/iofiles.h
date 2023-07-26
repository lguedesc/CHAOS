#ifndef IOFILES_H
#define IOFILES_H

#include <stdio.h>
#include <stdbool.h>
#include "defines.h"

// Handles the reading, creation and checking of files and directories
char *convert_dir(const char *rawdir);
void OSmkdir(char *path);
void create_dir(const char *pathname);
void directory_exists(const char *pathname);
bool file_exists(const char* filename);
FILE *create_output_file(char *name, const char *ext, const char *dir);
FILE *name_and_create_output_files(const char *systemname, const char *directory, const char *module, const char *ext);
char *get_input_filename(void);

// Write Results in file
void write_results(FILE *output_file, int dim, double t, double *x, int mode);
void write_results_lyap(FILE *output_file, int dim, double t, double *lambda, double *s_lambda, int mode);
void p_write_epbasin_results(FILE *output_file, double **results, int pixels, int dim);

// Write results in file for Harmonic Oscillator (HOS) Toolbox
void HOS_write_timeseries_results(FILE *output_file, int dim, double t, double *x, int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, ang_info *angles, int mode);
void HOS_write_poinc_results(FILE *output_file, int dim, double t, double *x, ang_info *angles, int mode);
void HOS_write_ftimeseries_results(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, ang_info *angles, int mode);
void HOS_write_bifurc_results(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int nrms, int *rmsindex, double *xrms, double *overallxrms,
                             int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, ang_info *angles, int mode);
void HOS_write_fbifurc_results(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *overallxmin, double *overallxmax, double *LE, int attractor, size_t npoinc, double **poinc, int nrms, int *rmsindex, double *xrms, double *overallxrms,
                              int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, ang_info *angles, int mode);
void HOS_p_write_dyndiag_results(FILE *output_file, int dim, int nrms, int *rmsindex, double **results, int pixels, int ncustomvalues, char **customnames, int nprintf, int *printfindex, ang_info *angles);
void HOS_p_write_fdyndiag_results(FILE *output_file, int dim, int nrms, int *rmsindex, double **results, int pixels, int ncustomvalues, char **customnames, int nprintf, int *printfindex, ang_info *angles);


#endif
