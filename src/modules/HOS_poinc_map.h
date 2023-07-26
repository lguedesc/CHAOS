#ifndef HOS_POINC_MAP_H
#define HOS_POINC_MAP_H

void HOS_poincaremap(char *funcname, unsigned int DIM, unsigned int nPar, ang_info *angles, char* outputname, void (*edosys)(int, double *, double, double *, double *));

#endif
