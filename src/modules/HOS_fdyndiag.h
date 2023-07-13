#ifndef HOS_FDYNDIAG_H
#define HOS_FDYNDIAG_H

void HOS_fdyndiag(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *), void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int));

#endif