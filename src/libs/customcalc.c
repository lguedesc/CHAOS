#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "energyharvest.h"

// Methods
void assign_names(char **strings, const int nvalues, char **names, size_t maxstrlen) {
    // Get every row of strings and copy to names
    for (int i = 0; i < nvalues; i++) {
        if (strlen(strings[i]) + 1 >= maxstrlen) {
            printf("  CUSTOM CALCULATIONS FUNCTION WARNING: One or more given names of the custom values are too big, please assign smaller names of maximum length of %zu before trying to run the program.\n", maxstrlen);
            printf("  Exiting Program...\n");
            exit(1);
        }
        strcpy(names[i], strings[i]);
    }
}

// Custom Calculations
void customcalc(double *x, double *par, double t, double *xrms, int N, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {
    return;
}

void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, int N, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {
    // Check if mode is equal to "table"
    if (mode == 0) {
        // Names to be printed in the output file
        char *names[] = {   "xb",
                            "dxb",
                            "ddxb",
                            "Sum(xb^2)",
                            "Sum(dxb^2)",
                            "Sum(ddxb^2)",
                            "TRdispl",
                            "TRvel",
                            "xbRMS",
                            "dxbRMS",
                            "ddxbRMS",
                            "PoutAvg",
                            "PinAvg",
                            "EffAvg"
                        };
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames, maxstrlen);
    }
    else if (mode == 1) {
        // Input Base Excitation Displacement
        customvalue[0] = par[1]*sin(par[0]*t);
        // Input Base Excitation Velocity
        customvalue[1] = par[1]*par[0]*cos(par[0]*t);
        // Input Base Excitation Acceleration
        customvalue[2] = par[1]*par[0]*par[0]*sin(par[0]*t);
        // Accumulate the value of the square of the Input Base Excitation Displacement
        customvalue[3] = RMS(&customvalue[3], customvalue[0], N, 0);
        // Accumulate the value of the square of the Input Base Excitation Velocity
        customvalue[4] = RMS(&customvalue[4], customvalue[1], N, 0);
        // Accumulate the value of the square of the Input Base Excitation Acceleration
        customvalue[5] = RMS(&customvalue[5], customvalue[2], N, 0);
        // Transmissibility of displacement
        customvalue[6] = x[0]/customvalue[0];
        // Transmissibility of velocity
        customvalue[7] = x[1]/customvalue[1];
    }
    // Check if mode is equal to "end"
    else if (mode == 2) {
        // RMS of the Input Base Excitation Displacement
        customvalue[8] = RMS(&customvalue[3], customvalue[0], N, 1);
        // RMS of the Input Base Excitation Velocity
        customvalue[9] = RMS(&customvalue[4], customvalue[1], N, 1);
        // RMS of the Input Base Excitation Acceleration
        customvalue[10] = RMS(&customvalue[5], customvalue[2], N, 1);
        // Peak to Peak Average Electrical Output Power
        customvalue[11] = (par[5]*par[6]/par[7])*xrms[2];
        // Peak to Peak Average Mechanical Input Power of the relative motion with respect to the base
        customvalue[12] = customvalue[10]*xrms[1];
        // Efficiency of Conversion n = pe/pm
        customvalue[13] = customvalue[11]/customvalue[12];
    }
    else {
        printf("DEBUG WARNING: Custom Function using mode = %d, please use 0, 1 or 2", mode);
    }
}