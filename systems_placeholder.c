#include <stdio.h>

void linear_osc(int dim, double *x, double t, double *par, double *f) {
    printf("Running linear_osc with dim=%d\n", dim);
}

void duffing(int dim, double *x, double t, double *par, double *f) {
    printf("Running duffing with dim=%d\n", dim);
}

void vanderpol(int dim, double *x, double t, double *par, double *f) {
    printf("Running vanderpol with dim=%d\n", dim);
}

void bistable_EH(int dim, double *x, double t, double *par, double *f) {
    printf("Running bistable_EH with dim=%d\n", dim);
}

void num_integrator(void (*odesys)(int, double *, double, double *, double *), double dim, double *x, double t, double *par, double *f) {
    odesys(dim, x, t, par, f);
}

void customcalc_duffing(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
    printf("Running customcalc_duffing_EH.\n");
}

void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
    printf("Running customcalc_bistable_EH\n");
}