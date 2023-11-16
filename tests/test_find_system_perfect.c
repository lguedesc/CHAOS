#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Define your ODE system functions
void duffing(double a, double b) {
    printf("Running duffing with a=%f, b=%f...\n", a, b);
}

void vanderpol(double a, double b) {
    printf("Running vanderpol with a=%f, b=%f...\n", a, b);
}

void lorenz(double a, double b) {
    printf("Running lorenz with a=%f, b=%f...\n", a, b);
}

// Define a struct to hold the system name and function pointer
typedef struct {
    const char* name;
    void (*func)(double, double);
} System;

// Define an array of systems
System systems[] = {
    {"duffing", duffing},
    {"vanderpol", vanderpol},
    {"lorenz", lorenz}
};

// Function to find a system by name and return the function pointer
void (*findSystem(const char *name))(double, double) {
    for (int i = 0; i < sizeof(systems)/sizeof(System); i++) {
        if (strcmp(name, systems[i].name) == 0) {
            return systems[i].func;
        }
    }
    return NULL;
}

// Your num_integrator function
void num_integrator(void (*odesys)(double, double), double a, double b) {
    odesys(a, b);
}

int main() {
    char *systemName = "duffing";
    void (*systemFunc)(double, double) = findSystem(systemName);
    if (systemFunc != NULL) {
        num_integrator(systemFunc, 1.0, 2.0);
    } else {
        printf("System not found: %s\n", systemName);
    }

    return 0;
}
