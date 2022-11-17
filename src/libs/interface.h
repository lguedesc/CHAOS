#include <stdio.h>
#include <stdlib.h>

// Main Interface
void clear_screen();
void partition(int mode, size_t maxlength);
void welcome_header(size_t maxlength);
void option_header(char *tbnames, size_t maxlength);
void mixed_header(char *names1, char *names2, size_t maxlength);
unsigned int choose_option(char **options, size_t n, size_t maxlength, char *category);
void invalid_option(unsigned int option, char* category, size_t maxlength);
void end_of_execution(size_t maxlength);
void identify_simulation(unsigned int toolbox, unsigned int *system, unsigned int *module, char **toolboxesNames, char **systemNames, char** moduleNames, size_t numofsystems, size_t maxlength, size_t numofmodules);
int int_length(int value);
// Simulation Prints
void fpartition(FILE *output_file, int mode, size_t maxlength);
void fwrite_module_and_system(FILE *output_file, char *funcname, char *modulename, size_t maxlength);

void write_initial_conditions(int dim, double *x, double t, size_t maxlength, double percname);
void fwrite_initial_conditions(FILE *output_file, int dim, double *x, double t, size_t maxlength, double percname);

void write_sys_parameters(int npar, double *par, size_t maxlength, double percname);
void fwrite_sys_parameters(FILE *output_file, int npar, double *par, size_t maxlength, double percname);

void write_RMS_calculations_info(int n, int *index, size_t maxlength, double percname);
void fwrite_RMS_calculations_info(FILE *output_file, int n, int *index, size_t maxlength, double percname);

void write_custom_info_calculations(int n, int nf, int *findex, int nscr, int *scrindex, size_t maxlength, double percname);
void fwrite_custom_info_calculations(FILE* output_file, int n, int nf, int *findex, int nscr, int *scrindex, size_t maxlength, double percname);

void write_prog_parameters_timeseries(int dim, int npar, int np, int ndiv, int trans, double h, size_t maxlength, double percname);
void fwrite_prog_parameters_timeseries(FILE* output_file, char *funcname, int dim, int npar, int np, int ndiv, int trans, double h, size_t maxlength, double percname);

void write_prog_parameters_ftimeseries(int dim, int npar, int maxper, int np, int ndiv, int trans, double h, size_t maxlength, double percname);
void fwrite_prog_parameters_ftimeseries(FILE* output_file, char *funcname, int dim, int npar, int maxper, int np, int ndiv, int trans, double h, size_t maxlength, double percname);

void write_prog_parameters_bifurcation(int dim, int npar, int np, int ndiv, int trans, size_t maxlength, double percname);
void fwrite_prog_parameters_bifurcation(FILE* output_file, char *funcname, int dim, int npar, int np, int ndiv, int trans, size_t maxlength, double percname);

void write_prog_parameters_fbifurcation(int dim, int npar, int np, int ndiv, int trans, int maxper, size_t maxlength, double percname);
void fwrite_prog_parameters_fbifurcation(FILE *output_file, char *funcname, int dim, int npar, int np, int ndiv, int trans, int maxper, size_t maxlength, double percname);

void write_prog_parameters_dyndiag(int dim, int npar, int np, int ndiv, int trans, size_t maxlength, double percname);
void fwrite_prog_parameters_dyndiag(FILE *output_file, char* funcname, int dim, int npar, int np, int ndiv, int trans, size_t maxlength, double percname);

void write_prog_parameters_fdyndiag(int dim, int npar, int np, int ndiv, int maxper, int trans, size_t maxlength, double percname);
void fwrite_prog_parameters_fdyndiag(FILE *output_file, char *funcname, int dim, int npar, int np, int ndiv, int maxper, int trans, size_t maxlength, double percname);

void write_bifurcation_info(double *parrange, int parindex, int bifmode, size_t maxlength, double percname);
void fwrite_bifurcation_info(FILE* output_file, double *parrange, int parindex, int bifmode, size_t maxlength, double percname);

void write_dyndiag_info(double *parrange, int indexX, int indexY, int bifmode, size_t maxlength, double percname);
void fwrite_dyndiag_info(FILE *output_file, double *parrange, int indexX, int indexY, int bifmode, size_t maxlength, double percname);

void print_RMS(int nRMS, int *rmsindex, double *xRMS, double *overallxRMS, size_t maxlength, double percname);
void fprint_RMS(FILE *output_file, int nRMS, int *rmsindex, double *xRMS, double *overallxRMS, size_t maxlength, double percname);

void print_customcalc(int nprintscr, int *printscrindex, double *customvalue, char **customnames, size_t maxlength, double percname);
void fprint_customcalc(FILE *output_file, int nprintscr, int *printscrindex, double *customvalue, char **customnames, size_t maxlength, double percname);

void print_minmax(double *xmin, double *xmax, double *overallxmin, double *overallxmax, int dim, size_t maxlength, double percname);
void fprint_minmax(FILE *output_file, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int dim, size_t maxlength, double percname);

void print_attractor(int attrac, int maxper, size_t maxlength, double percname);
void fprint_attractor(FILE *output_file, int attrac, int maxper, size_t maxlength, double percname);