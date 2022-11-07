#include <stdio.h>

void EH_timeseries(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, int, int, char **, size_t, double *, int));