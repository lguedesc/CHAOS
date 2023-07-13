#ifndef HOS_FFORCEDBASIN_H
#define HOS_FFORCEDBASIN_H

void HOS_fforcedbasin(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int));

#endif