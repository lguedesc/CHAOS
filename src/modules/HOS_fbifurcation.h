#ifndef HOS_FBIFURCATION_H
#define HOS_FBIFURCATION_H

void HOS_fbifurcation(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int));

#endif