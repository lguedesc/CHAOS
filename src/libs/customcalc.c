#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nlosc.h"
#include "nldyn.h"

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
void customcalc(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {
    return;
}

void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {
    // Mode to define names to be printed in the output file
    if (mode == 0) {    
        char *names[] = {   "ddx[0]",
                            "ddx[0]MIN",
                            "ddx[0]MAX",
                            "Sum(ddx[0]^2)",
                            "OVRLLddx[0]",
                            "xb",
                            "dxb",
                            "ddxb",
                            "OVRLLddx[0]MIN",
                            "OVRLLddx[0]MAX",
                            "Sum(OVRLLddx[0]^2)",
                            "Sum(OVRLLxb^2)",
                            "Sum(OVRLLdxb^2)",
                            "Sum(OVRLLddxb^2)",
                            "ddx[0]RMS",
                            "OVRLLddx[0]RMS",
                            "xbRMS",
                            "dxbRMS",
                            "ddxbRMS",
                            "PoutAvg"
                        };
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames, maxstrlen);
    }
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) { 
        // Acceleration of the system in steady state regime
        customvalue[0] = par[1]*par[0]*par[0]*sin(par[0]*t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] + par[5]*x[2];
        // Minimum value of the acceleration in steady state regime
        min_value(customvalue[0], &customvalue[1]);
        // Maximum Value of the Acceleration in steady state regime
        max_value(customvalue[0], &customvalue[2]);
        // Accumulate the value of the square of the Acceleration of the system in steady state regime
        customvalue[3] = RMS(&customvalue[3], customvalue[0], (int)(N*steadystateperc), 0);
    }
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        double ddxb = par[1]*par[0]*par[0]*sin(par[0]*t);
        // Overall Acceleration of the system 
        customvalue[4] = ddxb - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] + par[5]*x[2];
        // Overall Input Base Excitation Displacement
        customvalue[5] = par[1]*sin(par[0]*t);
        // Overall Input Base Excitation Velocity
        customvalue[6] = par[1]*par[0]*cos(par[0]*t);
        // Overall Input Base Excitation Acceleration
        customvalue[7] = -ddxb;
        // Overall minimum value of the acceleration of the system
        min_value(customvalue[4], &customvalue[8]);
        // Overall maximum value of the acceleration of the system
        max_value(customvalue[4], &customvalue[9]);
        // Accumulate the value of the square of the overall acceleration of the system
        customvalue[10] = RMS(&customvalue[10], customvalue[4], N, 0);
        // Accumulate the value of the square of the overall input base excitation displacement of the system
        customvalue[11] = RMS(&customvalue[11], customvalue[5], N, 0);
        // Accumulate the value of the square of the overall input base excitation velocity of the system
        customvalue[12] = RMS(&customvalue[12], customvalue[6], N, 0);
        // Accumulate the value of the square of the overall input base excitation acceleration of the system
        customvalue[13] = RMS(&customvalue[13], customvalue[7], N, 0);
    } 
    // Mode to perform calculations at the end of the time series
    else if (mode == 3) {
        // RMS acceleration of the system in the steady state regime
        customvalue[14] = RMS(&customvalue[3], customvalue[0], (int)(N*steadystateperc), 1);
        // Overall RMS acceleration of the system 
        customvalue[15] = RMS(&customvalue[10], customvalue[4], N, 1);
        // Overall RMS of the input base excitation displacement
        customvalue[16] = RMS(&customvalue[11], customvalue[5], N, 1);
        // Overall RMS of the input base excitation velocity
        customvalue[17] = RMS(&customvalue[12], customvalue[6], N, 1); 
        // Overall RMS of the input base excitation acceleration
        customvalue[18] = RMS(&customvalue[13], customvalue[7], N, 1);
        // Peak to Peak Average Electrical Output Power
        customvalue[19] = (par[5]*par[6]/par[7])*xrms[2];
    }
    else {
        printf("DEBUG WARNING: Custom Function using mode = %d, please use 0, 1, 2 or 3", mode);
    }
}



void customcalc_bistable_EH_old2(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {
    // Check if mode is equal to "table"
    if (mode == 0) {    // Mode to define names to be printed in the output file
        char *names[] = {   "ddx[0]",
                            "xb",
                            "dxb",
                            "ddxb",
                            "PinInst",
                            "PinInst(Sum)",
                            "PoutAvg",
                            "PinAvg",
                            "EffAvg",
                        };
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames, maxstrlen);
    }
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) { 
        // Acceleration of the system
        customvalue[0] = par[1]*par[0]*par[0]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] + par[5]*x[2];
        // Input Base Excitation Displacement
        customvalue[1] = par[1]*sin(par[0]*t);
        // Input Base Excitation Velocity
        customvalue[2] = par[1]*par[0]*cos(par[0]*t);
        // Input Base Excitation Acceleration
        customvalue[3] = -par[1]*par[0]*par[0]*sin(par[0]*t);
        // Instantaneous Input Power
        customvalue[4] = (customvalue[0] + customvalue[3])*customvalue[2];
        // Cumulative sum of Instantaneous Input Power
        customvalue[5] = customvalue[5] + customvalue[4];
    }
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        // Acceleration of the system over the entire time series
        customvalue[0] = par[1]*par[0]*par[0]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] + par[5]*x[2];
        // Input Base Excitation Displacement over the entire time series
        customvalue[1] = par[1]*sin(par[0]*t);
        // Input Base Excitation Velocity over the entire time series
        customvalue[2] = par[1]*par[0]*cos(par[0]*t); 
        // Input Base Excitation Acceleration over the entire time series
        customvalue[3] = -par[1]*par[0]*par[0]*sin(par[0]*t);
        // Instantaneous Input Power over the entire time series
        customvalue[4] = (customvalue[0] + customvalue[3])*customvalue[2];
        // Cumulative sum of Instantaneous Input Power over the entire time series
        customvalue[5] = customvalue[5] + customvalue[4];
    } 
    // Mode to perform calculations at the end of the time series
    else if (mode == 3) {
        // Peak to Peak Average Electrical Output Power
        customvalue[6] = (par[5]*par[6]/par[7])*xrms[2];
        // Peak to Peak Average Mechanical Input Power of the relative motion with respect to the base
        customvalue[7] = customvalue[5]/N;
        // Efficiency of Conversion n = pe/pm
        customvalue[8] = customvalue[6]/customvalue[7];
    }
    else {
        printf("DEBUG WARNING: Custom Function using mode = %d, please use 0, 1 or 2", mode);
    }
}

void customcalc_bistable_EH_old(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {
    // Check if mode is equal to "table"
    if (mode == 0) {
        // Names to be printed in the output file
        char *names[] = {   "xb",
                            "dxb",
                            "ddxb",
                            "Sum(xb^2)",
                            "Sum(dxb^2)",
                            "Sum(ddxb^2)",
                            "xbRMS",
                            "dxbRMS",
                            "ddxbRMS",
                            "PoutAvg",
                            "PinAvg",
                            "EffAvg",
                            "TRdispl(withRMS)",
                            "TRvel(withRMS)",
                            "TRdispl(PtoP)",
                            "TRvel(PtoP)",
                        };
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames, maxstrlen);
    }
    else if (mode == 1) {
        // Input Base Excitation Displacement
        customvalue[0] = (par[1]/par[0]*par[0])*sin(par[0]*t);
        // Input Base Excitation Velocity
        customvalue[1] = (par[1]/par[0])*cos(par[0]*t);
        // Input Base Excitation Acceleration
        customvalue[2] = -par[1]*sin(par[0]*t);
        // Accumulate the value of the square of the Input Base Excitation Displacement
        customvalue[3] = RMS(&customvalue[3], customvalue[0], N, 0);
        // Accumulate the value of the square of the Input Base Excitation Velocity
        customvalue[4] = RMS(&customvalue[4], customvalue[1], N, 0);
        // Accumulate the value of the square of the Input Base Excitation Acceleration
        customvalue[5] = RMS(&customvalue[5], customvalue[2], N, 0);
    }
    // Check if mode is equal to "end"
    else if (mode == 2) {
        // RMS of the Input Base Excitation Displacement
        customvalue[6] = RMS(&customvalue[3], customvalue[0], N, 1);
        // RMS of the Input Base Excitation Velocity
        customvalue[7] = RMS(&customvalue[4], customvalue[1], N, 1);
        // RMS of the Input Base Excitation Acceleration
        customvalue[8] = RMS(&customvalue[5], customvalue[2], N, 1);
        // Peak to Peak Average Electrical Output Power
        customvalue[9] = (par[5]*par[6]/par[7])*xrms[2];
        // Peak to Peak Average Mechanical Input Power of the relative motion with respect to the base
        customvalue[10] = customvalue[8]*xrms[1];
        // Efficiency of Conversion n = pe/pm
        customvalue[11] = customvalue[9]/customvalue[10];
        // Avg Transmissibility of displacement (with RMS)
        customvalue[12] = xrms[0]/customvalue[6];
        // Avg Transmissibility of velocity (with RMS)
        customvalue[13] = xrms[1]/customvalue[7];
        // Avg Transmissibility of displacement (Peak to Peak)
        customvalue[14] = (xmax[0] - xmin[0])/(par[1] - (-par[1]));
        // Avg Transmissibility of velocity (Peak to Peak)
        customvalue[15] = (xmax[1] - xmin[1])/(par[0]*(par[1] - (-par[1])));
    }
    else {
        printf("DEBUG WARNING: Custom Function using mode = %d, please use 0, 1 or 2", mode);
    }
}