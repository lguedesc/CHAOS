#ifndef SYSTEMS_PLACEHOLDER_H
#define SYSTEMS_PLACEHOLDER_H

void duffing(int dim, double *x, double t, double *par, double *f);
void bistable_EH(int dim, double *x, double t, double *par, double *f);

void num_integrator(void (*odesys)(int, double *, double, double *, double *), double dim, double *x, double t, double *par, double *f);

void customcalc_duffing(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode);

void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode);

#endif
