#include <stdio.h>

void duffing(double a, double b) {
    printf("Running duffing with a=%f, b=%f...\n", a, b);
}

void vanderpol(double a, double b) {
    printf("Running vanderpol with a=%f, b=%f...\n", a, b);
}

void lorenz(double a, double b) {
    printf("Running lorenz with a=%f, b=%f...\n", a, b);
}

void bistable_EH(double a, double b) {
    printf("Running bistable_EH with a=%f, b=%f...\n", a, b);
}

void num_integrator(void (*odesys)(double, double), double a, double b) {
    odesys(a, b);
}