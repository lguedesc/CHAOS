#include <stdio.h>

void OS_fdyndiag(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, int, double, int, char **, size_t, double *, int));