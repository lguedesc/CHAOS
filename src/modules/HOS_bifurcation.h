#ifndef HOS_BIFURCATION_H
#define HOS_BIFURCATION_H

void HOS_bifurcation(char *funcname, unsigned int DIM, unsigned int nPar, ang_info *angles, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int));

#endif