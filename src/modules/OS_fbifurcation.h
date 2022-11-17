#include <stdio.h>

void OS_fbifurcation(char *funcname, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, int, int, char **, size_t, double *, int));