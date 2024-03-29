#ifndef ODESYSTEMS_H
#define ODESYSTEMS_H

// General Nonlinear Systems
void lorenz(int dim, double *x, double t, double *par, double *f);
void lotka_volterra_predator_prey(int dim, double *x, double t, double *par, double *f);
void halvorsen(int dim, double *x, double t, double *par, double *f);

// Nonlinear Oscillators
void falksma(int dim, double *x, double t, double *par, double *f);
void vanderpol(int dim, double *x, double t, double *par, double *f);
void pendulum(int dim, double *x, double t, double *par, double *f);
void duffing(int dim, double *x, double t, double *par, double *f);
void linear_oscillator(int dim, double *x, double t, double *par, double *f);
void duffing_2DoF(int dim, double *x, double t, double *par, double *f);
void duffing_vanderpol(int dim, double *x, double t, double *par, double *f);
void linear_oscillator_2DoF(int dim, double *x, double t, double *par, double *f);

// Mechanical Energy Harvesters
void bistable_EH(int dim, double *x, double t, double *par, double *f);
void tristable_EH(int dim, double *x, double t, double *par, double *f);
void tetrastable_EH(int dim, double *x, double t, double *p, double *f);
void pend_oscillator_EH(int dim, double *x, double t, double *par, double *f);
void pend_oscillator_wout_pend_EH(int dim, double *x, double t, double *par, double *f);
void duffing_2DoF_EH(int dim, double *x, double t, double *par, double *f);
void linear_2DoF_EH(int dim, double *x, double t, double *par, double *f);
void linear_EMEH(int dim, double *x, double t, double *par, double *f);
void pendulum_EMEH(int dim, double *x, double t, double *p, double *f);
void linear_oscillator_gravity(int dim, double *x, double t, double *p, double *f);
void multidirectional_hybrid_EH(int dim, double *x, double t, double *p, double *f);
void multidirectional_hybrid_EH_zero_pend_length(int dim, double *x, double t, double *p, double *f);
void pendulum_EMEH_dimensional(int dim, double *x, double t, double *p, double *f);
void multidirectional_hybrid_EH_coupling_ratio(int dim, double *x, double t, double *p, double *f);

// Adeodato SMA Model
void adeodato_sma_oscillator(int dim, double *x, double t, double *par, double *f);

// Not Implemented
void duffing_cldyn(int dim, double *x, double t, double *par, double *f);

// Automatically Added
void chuas_circuit(int dim, double *x, double t, double *par, double *f);

#endif
