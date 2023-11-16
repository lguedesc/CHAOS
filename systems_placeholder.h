#ifndef SYSTEMS_PLACEHOLDER_H
#define SYSTEMS_PLACEHOLDER_H

void duffing(double a, double b);
void vanderpol(double a, double b);
void lorenz(double a, double b);
void bistable_EH(double a, double b);

void num_integrator(void (*odesys)(double, double), double a, double b);

#endif
