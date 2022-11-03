#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "energyharvest.h"

void customcalc(double *x, double *par, double t, double *xrms, int N, int ncustomvalues, char *customnames[], double *customvalue, char* mode) {
    return;
}

void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, int N, int ncustomvalues, char *customnames[], double *customvalue, char* mode) {
    // Check if mode is equal to "table"
    if (strcmp(mode, "table") == 0) {
        // Input Base Excitation Displacement
        customvalue[0] = par[1]*sin(par[0]*t);
        // Input Base Excitation Velocity
        customvalue[1] = par[1]*par[0]*cos(par[0]*t);
        // Input Base Excitation Acceleration
        customvalue[2] = par[1]*par[0]*par[0]*sin(par[0]*t);
        // Accumulate the value of the square of the Input Base Excitation Displacement
        customvalue[3] = customvalue[3] + (customvalue[0]*customvalue[0]);
        // Accumulate the value of the square of the Input Base Excitation Velocity
        customvalue[4] = customvalue[4] + (customvalue[1]*customvalue[1]);
        // Accumulate the value of the square of the Input Base Excitation Acceleration
        customvalue[5] = customvalue[5] + (customvalue[2]*customvalue[2]);
    }
    // Check if mode is equal to "end"
    else if (strcmp(mode, "end") == 0) {
        // RMS of the Input Base Excitation Displacement
        customvalue[6] = sqrt(customvalue[3] / N);
        // RMS of the Input Base Excitation Velocity
        customvalue[7] = sqrt(customvalue[4] / N);
        // RMS of the Input Base Excitation Acceleration
        customvalue[8] = sqrt(customvalue[5] / N);
        // Peak to Peak Average Electrical Output Power
        customvalue[9] = (par[5]*par[6]/par[7])*xrms[2];
        // Peak to Peak Average Mechanical Input Power of the relative motion with respect to the base
        customvalue[10] = customvalue[8]*xrms[1];
        // Efficiency of Conversion n = pe/pm
        customvalue[11] = customvalue[9]/customvalue[10];
    }
    // Check if mode is equal to "names"
    else if (strcmp(mode, "names") == 0) {
        // Names to be printed in the output file
        customnames[0] = "xb";
        customnames[1] = "dxb";
        customnames[2] = "ddxb";
        customnames[3] = "Sum(xb^2)";
        customnames[4] = "Sum(dxb^2)";
        customnames[5] = "Sum(ddxb^2)";
        customnames[6] = "xbRMS";
        customnames[7] = "dxbRMS";
        customnames[8] = "ddxbRMS";
        customnames[9] = "Pout";
        customnames[10] = "Pin";
        customnames[11] = "Eff";
    }
}