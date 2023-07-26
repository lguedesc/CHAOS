#ifndef ODESOLVERS_H
#define ODESOLVERS_H

void rk4(int dim, double *x, double t, double h, double *par, double *f, void (*edosys)(int, double *, double, double *, double *));

#endif
