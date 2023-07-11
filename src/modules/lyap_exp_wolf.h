#ifndef LYAP_EXP_WOLF_H
#define LYAP_EXP_WOLF_H

#include <stdio.h>

void lyapunov_exp_wolf(char *funcname, unsigned int DIM, unsigned int nPar, char* outputname, void (*edosys)(int, double *, double, double *, double *));

#endif