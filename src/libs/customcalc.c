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

void spins(double previous_angle, double current_angle, double *positive_spin, double *negative_spin) {
    if (current_angle > previous_angle) {
        (*positive_spin) = (*positive_spin) + fabs(current_angle);
    }
    else if (current_angle < previous_angle) {
        (*negative_spin) = (*negative_spin) + fabs(current_angle);
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

void customcalc_pend_oscillator_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {
    /* OMEGA   = par[0]   |   zeta_z    = par[5]   |   l         = par[10]   |   chi_PZ = par[15]       |   x[0] = x       |   x[5] = dphi/dt
       gamma   = par[1]   |   zeta_t    = par[6]   |   varphi_PZ = par[11]   |   chi_EM = par[16]       |   x[1] = dx/dt   |   x[6] = v
       mu      = par[2]   |   OMEGA_s   = par[7]   |   kappa_PZ  = par[12]   |                          |   x[2] = z       |   x[7] = i
       rho     = par[3]   |   OMEGA_phi = par[8]   |   varphi_EM = par[13]   |                          |   x[3] = dz/dt   |
       zeta_x  = par[4]   |   OMEGA_t   = par[9]   |   kappa_EM  = par[14]   |                          |   x[4] = phi     |                   */
    
    const double pi = 4 * atan(1);
    double rb = par[1]*sin(par[0]*t);
    double drb = par[1]*par[0]*cos(par[0]*t);
    double ddrb = -par[1]*par[0]*par[0]*sin(par[0]*t);
    // Mode to define names to be printed in the output file
    if (mode == 0) {    
        char *names[] = {   "ddX", "ddZ", "ddPhi", "Xcm", "Zcm", "dXcm", "dZcm", "ddXcm", "ddZcm",
                            "ddX_MIN", "ddZ_MIN", "ddPhi_MIN",
                            "ddX_MAX", "ddZ_MAX", "ddPhi_MAX",
                            "Xcm_MIN", "Zcm_MIN", "dXcm_MIN", "dZcm_MIN", "ddXcm_MIN", "ddZcm_MIN",
                            "Xcm_MAX", "Zcm_MAX", "dXcm_MAX", "dZcm_MAX", "ddXcm_MAX", "ddZcm_MAX",
                            "Sum(ddX^2)", "Sum(ddZ^2)", "Sum(ddPhi^2)", 
                            "Sum(Xcm^2)", "Sum(Zcm^2)", "Sum(dXcm^2)", "Sum(dZcm^2)", "Sum(ddXcm^2)", "Sum(ddZcm^2)",
                            "OVRLL_ddX", "OVRLL_ddZ", "OVRLLddPhi", "OVRLL_Xcm", "OVRLL_Zcm", "OVRLL_dXcm", "OVRLL_dZcm", "OVRLL_ddXcm", "OVRLL_ddZcm", "OVRLL_Xb", "OVRLL_Zb", "OVRLL_dXb", "OVRLL_dZb", "OVRLL_ddXb", "OVRLL_ddZb",
                            "OVRLL_ddX_MIN", "OVRLL_ddZ_MIN", "OVRLL_ddPhi_MIN", 
                            "OVRLL_ddX_MAX", "OVRLL_ddZ_MAX", "OVRLL_ddPhi_MAX",
                            "OVRLL_Xcm_MIN", "OVRLL_Zcm_MIN", "OVRLL_dXcm_MIN", "OVRLL_dZcm_MIN", "OVRLL_ddXcm_MIN", "OVRLL_ddZcm_MIN",
                            "OVRLL_Xcm_MAX", "OVRLL_Zcm_MAX", "OVRLL_dXcm_MAX", "OVRLL_dZcm_MAX", "OVRLL_ddXcm_MAX", "OVRLL_ddZcm_MAX",
                            "Sum(OVRLL_ddX^2)", "Sum(OVRLL_ddZ^2)", "Sum(OVRLL_ddPhi^2)", 
                            "Sum(OVRLL_Xcm^2)", "Sum(OVRLL_Zcm^2)", "Sum(OVRLL_dXcm^2)", "Sum(OVRLL_dZcm^2)", "Sum(OVRLL_ddXcm^2)", "Sum(OVRLL_ddZcm^2)",
                            "Sum(OVRLL_Xb)", "Sum(OVRLL_Zb)", "Sum(OVRLL_dXb)", "Sum(OVRLL_dZb)", "Sum(OVRLL_ddXb)", "Sum(OVRLL_ddZb)",  
                            "ddX_RMS", "ddZ_RMS", "ddPhi_RMS",
                            "Xcm_RMS", "Zcm_RMS", "dXcm_RMS", "dZcm_RMS", "ddXcm_RMS", "ddZcm_RMS",
                            "OVRLL_ddX_RMS", "OVRLL_ddZ_RMS", "OVRLL_ddPhi_RMS",
                            "OVRLL_Xcm_RMS", "OVRLL_Zcm_RMS", "OVRLL_dXcm_RMS", "OVRLL_dZcm_RMS", "OVRLL_ddXcm_RMS", "OVRLL_ddZcm_RMS",
                            "OVRLL_Xb_RMS", "OVRLL_Zb_RMS", "OVRLL_dXb_RMS", "OVRLL_dZb_RMS", "OVRLL_ddXb_RMS", "OVRLL_ddZb_RMS",
                            "PoutPZ_Avg", "PoutEM_Avg"};
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames, maxstrlen);
    }
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) { 
        // X Acceleration of the system in steady state regime ("ddX")
        customvalue[0] = (1/(1 + par[3]))*(-(1 + par[3]*cos(x[4])*cos(x[4]))*(2*par[4]*x[1] + par[7]*par[7]*x[0]) + (par[3]/2)*(2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*sin(x[4]))
                          + (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*cos(x[4]) -ddrb*sin((pi/180)*par[2]);
        // Z Acceleration of the system in steady state regime ("ddZ")
        customvalue[1] = (1/(1 + par[3]))*(-(1 + par[3]*sin(x[4])*sin(x[4]))*(2*par[5]*x[3] + x[2] - par[15]*x[6]) + (par[3]/2)*(2*par[4]*x[1] + par[7]*par[7]*x[0])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*cos(x[4]))
                          - (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*sin(x[4]) -ddrb*cos((pi/180)*par[2]);
        // Phi Acceleration of the system in steady state regime ("ddPhi")
        customvalue[2] = -(1 + par[3])*(2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7]) + (1/par[10])*((2*par[4]*x[1] + par[7]*par[7]*x[0])*cos(x[4]) - (2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(x[4]));
        // Displacement of the center of mass in X direction ("Xcm")
        customvalue[3] = (x[0] + par[3]*par[10]*sin(x[4]))/(1 + par[3]);
        // Displacement of the center of mass in Z direction ("Zcm")
        customvalue[4] = (x[2] + par[3]*par[10]*cos(x[4]))/(1 + par[3]);
        // Velocity of the center of mass in X direction ("dXcm")
        customvalue[5] = (x[1] + par[3]*par[10]*x[5]*cos(x[4]))/(1 + par[3]);
        // Velocity of the center of mass in Z direction ("dZcm")
        customvalue[6] = (x[3] - par[3]*par[10]*x[5]*sin(x[4]))/(1 + par[3]);
        // Acceleration of the center of mass in X direction ("ddXcm")
        customvalue[7] = (customvalue[0] + par[3]*par[10]* (customvalue[2]*cos(x[4]) - x[5]*x[5]*sin(x[4])))/(1 + par[3]);
        // Acceleration of the center of mass in Z direction ("ddZcm")
        customvalue[8] = (customvalue[1] - par[3]*par[10]* (customvalue[2]*sin(x[4]) + x[5]*x[5]*cos(x[4])))/(1 + par[3]);
        // Minimum values of the accelerations in steady state regime ("ddX_MIN", "ddZ_MIN", "ddPhi_MIN")
        for (int i = 0; i < 3; i++) { // From customvalue[9] to customvalue[11]
            min_value(customvalue[i], &customvalue[9+i]);
        }
        // Maximum values of the accelerations in steady state regime ("ddX_MAX", "ddZ_MAX", "ddPhi_MAX")
        for (int i = 0; i < 3; i++) { // From customvalue[12] to customvalue[14]
            max_value(customvalue[i], &customvalue[12+i]);
        }
        // Minimum values of the values of the center of mass in steady state regime ("Xcm_MIN", "Zcm_MIN", "dXcm_MIN", "dZcm_MIN", "ddXcm_MIN", "ddZcm_MIN")
        for (int i = 0; i < 6; i++) { // From customvalue[15] to customvalue[20]
            min_value(customvalue[3+i], &customvalue[15+i]);
        }
        // Maximum values of the values of the center of mass in steady state regime ("Xcm_MAX", "Zcm_MAX", "dXcm_MAX", "dZcm_MAX", "ddXcm_MAX", "ddZcm_MAX")
        for (int i = 0; i < 6; i++) { // From customvalue[21] to customvalue[26]
            max_value(customvalue[3+i], &customvalue[21+i]);
        }
        // Accumulate the value of the square of accelerations of the system in steady state regime ("Sum(ddX^2)", "Sum(ddZ^2)", "Sum(ddPhi^2)")
        for (int i = 0; i < 3; i++) { // From customvalue[27] to customvalue[29]
            customvalue[27+i] = RMS(&customvalue[27+i], customvalue[i], (int)(N*steadystateperc), 0);
        }
        // Accumulate the value of the square of values of the center of mass in steady state regime ("Sum(Xcm^2)", "Sum(Zcm^2)", "Sum(dXcm^2)", "Sum(dZcm^2)", "Sum(ddXcm^2)", "Sum(ddZcm^2)")
        for (int i = 0; i < 6; i++) { // From customvalue[30] to customvalue[35]
            customvalue[30+i] = RMS(&customvalue[30+i], customvalue[3+i], (int)(N*steadystateperc), 0);
        }     
    }
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        // Overall X Acceleration of the system ("OVRLL_ddX")
        customvalue[36] = (1/(1 + par[3]))*(-(1 + par[3]*cos(x[4])*cos(x[4]))*(2*par[4]*x[1] + par[7]*par[7]*x[0]) + (par[3]/2)*(2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*sin(x[4]))
                           + (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*cos(x[4]) -ddrb*sin((pi/180)*par[2]);
        // Overall Z Acceleration of the system ("OVRLL_ddZ")
        customvalue[37] = (1/(1 + par[3]))*(-(1 + par[3]*sin(x[4])*sin(x[4]))*(2*par[5]*x[3] + x[2] - par[15]*x[6]) + (par[3]/2)*(2*par[4]*x[1] + par[7]*par[7]*x[0])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*cos(x[4]))
                          - (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*sin(x[4]) -ddrb*cos((pi/180)*par[2]);
        // Overall Phi Acceleration of the system ("OVRLL_ddPhi")
        customvalue[38] = -(1 + par[3])*(2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7]) + (1/par[10])*((2*par[4]*x[1] + par[7]*par[7]*x[0])*cos(x[4]) - (2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(x[4]));
        // Overall displacement of the center of mass in X direction ("OVRLL_Xcm")
        customvalue[39] = (x[0] + par[3]*par[10]*sin(x[4]))/(1 + par[3]); 
        // Overall displacement of the center of mass in Z direction ("OVRLL_Zcm")
        customvalue[40] = (x[2] + par[3]*par[10]*cos(x[4]))/(1 + par[3]);
        // Overall velocity of the center of mass in X direction ("OVRLL_dXcm")
        customvalue[41] = (x[1] + par[3]*par[10]*x[5]*cos(x[4]))/(1 + par[3]);
        // Overall velocity of the center of mass in Z direction ("OVRLL_dZcm")
        customvalue[42] = (x[3] - par[3]*par[10]*x[5]*sin(x[4]))/(1 + par[3]);
        // Overall Acceleration of the center of mass in X direction ("OVRLL_ddXcm")
        customvalue[43] = (customvalue[36] + par[3]*par[10]*(customvalue[38]*cos(x[4]) - x[5]*x[5]*sin(x[4])))/(1 + par[3]);
        // Overall Acceleration of the center of mass in Z direction ("OVRLL_ddZcm")
        customvalue[44] = (customvalue[37] - par[3]*par[10]*(customvalue[38]*sin(x[4]) + x[5]*x[5]*cos(x[4])))/(1 + par[3]);
        // Overall Input Base Excitation Displacement in X direction ("OVRLL_Xb")
        customvalue[45] = rb*sin((pi/180)*par[2]);
        // Overall Input Base Excitation Displacement in Z direction ("OVRLL_Zb")
        customvalue[46] = rb*cos((pi/180)*par[2]);
        // Overall Input Base Excitation Velocity in X direction ("OVRLL_dXb")
        customvalue[47] = drb*sin((pi/180)*par[2]);
        // Overall Input Base Excitation Velocity in Z direction ("OVRLL_dZb")
        customvalue[48] = drb*cos((pi/180)*par[2]);
        // Overall Input Base Excitation Acceleration in X direction ("OVRLL_ddXb")
        customvalue[49] = ddrb*sin((pi/180)*par[2]);
        // Overall Input Base Excitation Acceleration in Z direction ("OVRLL_ddZb")
        customvalue[50] = ddrb*cos((pi/180)*par[2]);
        // Overall Minimum values of the accelerations of the system ("OVRLL_ddX_MIN", "OVRLL_ddZ_MIN", "OVRLL_ddPhi_MIN")
        for (int i = 0; i < 3; i++) {  // From customvalue[51] to customvalue[53]
            min_value(customvalue[36+i], &customvalue[51+i]);
        }
        // Overall Maximum values of the accelerations of the system ("OVRLL_ddX_MAX", "OVRLL_ddZ_MAX", "OVRLL_ddPhi_MAX")
        for (int i = 0; i < 3; i++) {  // From customvalue[54] to customvalue[56]
            max_value(customvalue[36+i], &customvalue[54+i]);
        }
        // Overall Minimum values of the values of the center of mass ("OVRLL_Xcm_MIN", "OVRLL_Zcm_MIN", "OVRLL_dXcm_MIN", "OVRLL_dZcm_MIN", "OVRLL_ddXcm_MIN", "OVRLL_ddZcm_MIN")
        for (int i = 0; i < 6; i++) { // From customvalue[57] to customvalue[62]
            min_value(customvalue[39+i], &customvalue[57+i]);
        }
        // Overall Maximum values of the values of the center of mass ("OVRLL_Xcm_MAX", "OVRLL_Zcm_MAX", "OVRLL_dXcm_MAX", "OVRLL_dZcm_MAX", "OVRLL_ddXcm_MAX", "OVRLL_ddZcm_MAX")
        for (int i = 0; i < 6; i++) { // From customvalue[63] to customvalue[68]
            min_value(customvalue[39+i], &customvalue[63+i]);
        }
        // Accumulate the value of the square of the overall accelerations of the system ("Sum(OVRLL_ddX^2)", "Sum(OVRLL_ddZ^2)", "Sum(OVRLL_ddPhi^2)")
        for (int i = 0; i < 3; i++) { // From customvalue[69] to customvalue[71]
            customvalue[69+i] = RMS(&customvalue[69+i], customvalue[36+i], N, 0);
        }
        // Accumulate the value of the square of the overall values of the center of mass ("Sum(OVRLL_Xcm^2)", "Sum(OVRLL_Zcm^2)", "Sum(OVRLL_dXcm^2)", "Sum(OVRLL_dZcm^2)", "Sum(OVRLL_ddXcm^2)", "Sum(OVRLL_ddZcm^2)")
        for (int i = 0; i < 6; i++) { // From customvalue[72] to customvalue[77]
            customvalue[72+i] = RMS(&customvalue[72+i], customvalue[39+i], N, 0);
        }
        // Accumulate the value of the square of the overall input base excitation values ("Sum(OVRLL_Xb)", "Sum(OVRLL_Zb)", "Sum(OVRLL_dXb)", "Sum(OVRLL_dZb)", "Sum(OVRLL_ddXb)", "Sum(OVRLL_ddZb)")
        for (int i = 0; i < 6; i++) { // From customvalue[78] to customvalue[83]
            customvalue[78+i] = RMS(&customvalue[78+i], customvalue[45+i], N, 0);
        }
    } 
    // Mode to perform calculations at the end of the time series    
    else if (mode == 3) {
        // RMS acceleration of the system in steady state regime ("ddX_RMS", "ddZ_RMS", "ddPhi_RMS")
        for (int i = 0; i < 3; i++) { // From customvalue[84] to customvalue[86]
            customvalue[84+i] = RMS(&customvalue[27+i], customvalue[i], (int)(N*steadystateperc), 1);
        }
        //RMS of the center of mass in steady state regime ("Xcm_RMS", "Zcm_RMS", "dXcm_RMS", "dZcm_RMS", "ddXcm_RMS", "ddZcm_RMS")
        for (int i = 0; i < 6; i++) { // From customvalue[87] to customvalue[92]
            customvalue[87+i] = RMS(&customvalue[30+i], customvalue[3+i], (int)(N*steadystateperc), 1);
        }     
        // Overall RMS acceleration of the system ("OVRLL_ddX_RMS", "OVRLL_ddZ_RMS", "OVRLL_ddPhi_RMS")
        for (int i = 0; i < 3; i++) { // From customvalue[93] to customvalue[95]
            customvalue[93+i] = RMS(&customvalue[69+i], customvalue[36+i], N, 1);
        }
        // Overall RMS of the center of mass ("OVRLL_Xcm_RMS", "OVRLL_Zcm_RMS", "OVRLL_dXcm_RMS", "OVRLL_dZcm_RMS", "OVRLL_ddXcm_RMS", "OVRLL_ddZcm_RMS")
        for (int i = 0; i < 6; i++) { // From customvalue[96] to customvalue[101]
            customvalue[96+i] = RMS(&customvalue[72+i], customvalue[39+i], N, 1);
        }
        // Overall RMS of the input base excitation ("OVRLL_Xb_RMS", "OVRLL_Zb_RMS", "OVRLL_dXb_RMS", "OVRLL_dZb_RMS", "OVRLL_ddXb_RMS", "OVRLL_ddZb_RMS")
        for (int i = 0; i < 6; i++) { // From customvalue[102] to customvalue[107]
            customvalue[102+i] = RMS(&customvalue[78+i], customvalue[45+i], N, 1);
        }
        // Peak to Peak Average Electrical Output Power of the Piezoelectric Element
        customvalue[108] = (par[15]*par[11]/par[12])*xrms[6];
        // Peak to Peak Average Electrical Output Power of the Electromagnetic Converter
        customvalue[109] = (par[16]*par[13]/par[14])*xrms[7];
        // !!!! TRY TO CALCULATE THE NUMBER OF FLIPS OF THE PENDULUM !!!!
        // !!!! ALSO DETERMINE THE AMOUNT OF ELAPSED TIME THE PENDULUM FLIPS OVER !!!!
    }
    else {
        printf("DEBUG WARNING: Custom Function using mode = %d, please use 0, 1, 2 or 3", mode);
    }
}



void customcalc_pend_oscillator_EH_old(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {
    /* OMEGA   = par[0]   |   zeta_z    = par[5]   |   l         = par[10]   |   chi_PZ = par[15]       |   x[0] = x       |   x[5] = dphi/dt
       gamma   = par[1]   |   zeta_t    = par[6]   |   varphi_PZ = par[11]   |   chi_EM = par[16]       |   x[1] = dx/dt   |   x[6] = v
       mu      = par[2]   |   OMEGA_s   = par[7]   |   kappa_PZ  = par[12]   |                          |   x[2] = z       |   x[7] = i
       rho     = par[3]   |   OMEGA_phi = par[8]   |   varphi_EM = par[13]   |                          |   x[3] = dz/dt   |
       zeta_x  = par[4]   |   OMEGA_t   = par[9]   |   kappa_EM  = par[14]   |                          |   x[4] = phi     |                   */
    
    const double pi = 4 * atan(1);
    double rb = par[1]*sin(par[0]*t);
    double drb = par[1]*par[0]*cos(par[0]*t);
    double ddrb = -par[1]*par[0]*par[0]*sin(par[0]*t);
    // Mode to define names to be printed in the output file
    if (mode == 0) {    
        char *names[] = {   "ddX", "ddZ", "ddPhi", 
                            "Xcm", "Zcm", 
                            "dXcm", "dZcm",
                            "ddXcm", "ddZcm",
                            "ddX_MIN", "ddX_MAX",
                            "ddZ_MIN", "ddZ_MAX",
                            "ddPhi_MIN", "ddPhi_MAX",
                            "Sum(ddX^2)", "Sum(ddZ^2)", "Sum(ddPhi^2)",
                            "Sum(Xcm^2)", "Sum(Zcm^2)",
                            "Sum(dXcm^2)","Sum(dZcm^2)",
                            "Sum(ddXcm^2)","Sum(ddZcm^2)",
                            "OVRLL_ddX", "OVRLL_ddZ", "OVRLL_ddPhi",
                            "OVRLL_Xcm", "OVRLL_Zcm",
                            "OVRLL_dXcm", "OVRLL_dZcm",
                            "OVRLL_ddXcm", "OVRLL_ddZcm",
                            "OVRLL_xb", "OVRLL_zb",
                            "OVRLL_dxb", "OVRLL_dzb",
                            "OVRLL_ddxb", "OVRLL_ddzb",    
                        };
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames, maxstrlen);
    }
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) { 
        // X Acceleration of the system in steady state regime
        customvalue[0] = (1/(1 + par[3]))*(-(1 + par[3]*cos(x[4])*cos(x[4]))*(2*par[4]*x[1] + par[7]*par[7]*x[0]) + (par[3]/2)*(2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*sin(x[4]))
                          + (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*cos(x[4]) -ddrb*sin((pi/180)*par[2]);
        // Z Acceleration of the system in steady state regime
        customvalue[1] = (1/(1 + par[3]))*(-(1 + par[3]*sin(x[4])*sin(x[4]))*(2*par[5]*x[3] + x[2] - par[15]*x[6]) + (par[3]/2)*(2*par[4]*x[1] + par[7]*par[7]*x[0])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*cos(x[4]))
                          - (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*sin(x[4]) -ddrb*cos((pi/180)*par[2]);
        // Phi Acceleration of the system in steady state regime
        customvalue[2] = -(1 + par[3])*(2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7]) + (1/par[10])*((2*par[4]*x[1] + par[7]*par[7]*x[0])*cos(x[4]) - (2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(x[4]));
        // Displacement of the center of mass in X direction
        customvalue[3] = (x[0] + par[3]*par[10]*sin(x[4]))/(1 + par[3]);
        // Displacement of the center of mass in Z direction
        customvalue[4] = (x[2] + par[3]*par[10]*cos(x[4]))/(1 + par[3]);
        // Velocity of the center of mass in X direction
        customvalue[5] = (x[1] + par[3]*par[10]*x[5]*cos(x[4]))/(1 + par[3]);
        // Velocity of the center of mass in Z direction
        customvalue[6] = (x[3] - par[3]*par[10]*x[5]*sin(x[4]))/(1 + par[3]);
        // Acceleration of the center of mass in X direction
        customvalue[7] = (customvalue[0] + par[3]*par[10]* (customvalue[2]*cos(x[4]) - x[5]*x[5]*sin(x[4])))/(1 + par[3]);
        // Acceleration of the center of mass in Z direction
        customvalue[8] = (customvalue[1] - par[3]*par[10]* (customvalue[2]*sin(x[4]) + x[5]*x[5]*cos(x[4])))/(1 + par[3]);
        // Minimum value of the X acceleration in steady state regime
        min_value(customvalue[0], &customvalue[9]);
        // Maximum Value of the X acceleration in steady state regime
        max_value(customvalue[0], &customvalue[10]);
        // Minimum value of the Z acceleration in steady state regime
        min_value(customvalue[1], &customvalue[11]);
        // Maximum value of the Z acceleration in steady state regime
        max_value(customvalue[1], &customvalue[12]);
        // Minimum value of the Phi acceleration in steady state regime
        min_value(customvalue[2], &customvalue[13]);
        // Maximum value of the Phi acceleration in steady state regime
        max_value(customvalue[2], &customvalue[14]);
        // Minimum value of the displacement of the center of mass in X direction in steady state regime
        // Maximum value of the displacement of the center of mass in Z direction in steady state regime
        // Minimum value of the velocity of the center of mass in X direction in steady state regime
        // Maximum value of the velocity of the center of mass in Z direction in steady state regime
        // Minimum value of the acceleration of the center of mass in X direction in steady state regime
        // Maximum value of the acceleration of the center of mass in Z direction in steady state regime
        // Accumulate the value of the square of the X Acceleration of the system in steady state regime
        customvalue[15] = RMS(&customvalue[15], customvalue[0], (int)(N*steadystateperc), 0);
        // Accumulate the value of the square of the Z Acceleration of the system in steady state regime
        customvalue[16] = RMS(&customvalue[16], customvalue[1], (int)(N*steadystateperc), 0);
        // Accumulate the value of the square of the Phi Acceleration of the system in steady state regime
        customvalue[17] = RMS(&customvalue[17], customvalue[2], (int)(N*steadystateperc), 0);
        // Accumulate the value of the square of the displacement of the center of mass in X direction in steady state regime
        customvalue[18] = RMS(&customvalue[18], customvalue[3], (int)(N*steadystateperc), 0);
        // Accumulate the value of the square of the displacement of the center of mass in Z direction in steady state regime
        customvalue[19] = RMS(&customvalue[19], customvalue[4], (int)(N*steadystateperc), 0);    
        // Accumulate the value of the square of the velocity of the center of mass in X direction in steady state regime            
        customvalue[20] = RMS(&customvalue[20], customvalue[5], (int)(N*steadystateperc), 0);
        // Accumulate the value of the square of the velocity of the center of mass in Z direction in steady state regime
        customvalue[21] = RMS(&customvalue[21], customvalue[6], (int)(N*steadystateperc), 0);
        // Accumulate the value of the square of the acceleration of the center of mass in X direction in steady state regime
        customvalue[22] = RMS(&customvalue[22], customvalue[7], (int)(N*steadystateperc), 0);
        // Accumulate the value of the square of the acceleration of the center of mass in Z direction in steady state regime
        customvalue[23] = RMS(&customvalue[23], customvalue[8], (int)(N*steadystateperc), 0);
    }
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        // Overall X Acceleration of the system 
        customvalue[24] = (1/(1 + par[3]))*(-(1 + par[3]*cos(x[4])*cos(x[4]))*(2*par[4]*x[1] + par[7]*par[7]*x[0]) + (par[3]/2)*(2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*sin(x[4]))
                           + (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*cos(x[4]) -ddrb*sin((pi/180)*par[2]);
        // Overall Z Acceleration of the system
        customvalue[25] = (1/(1 + par[3]))*(-(1 + par[3]*sin(x[4])*sin(x[4]))*(2*par[5]*x[3] + x[2] - par[15]*x[6]) + (par[3]/2)*(2*par[4]*x[1] + par[7]*par[7]*x[0])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*cos(x[4]))
                          - (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*sin(x[4]) -ddrb*cos((pi/180)*par[2]);
        // Overall Phi Acceleration of the system
        customvalue[26] = -(1 + par[3])*(2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7]) + (1/par[10])*((2*par[4]*x[1] + par[7]*par[7]*x[0])*cos(x[4]) - (2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(x[4]));
        // Overall displacement of the center of mass in X direction
        customvalue[27] = (x[0] + par[3]*par[10]*sin(x[4]))/(1 + par[3]);
        // Overall displacement of the center of mass in Z direction
        customvalue[28] = (x[2] + par[3]*par[10]*cos(x[4]))/(1 + par[3]);
        // Overall velocity of the center of mass in X direction
        customvalue[29] = (x[1] + par[3]*par[10]*x[5]*cos(x[4]))/(1 + par[3]);
        // Overall velocity of the center of mass in Z direction
        customvalue[30] = (x[3] - par[3]*par[10]*x[5]*sin(x[4]))/(1 + par[3]);
        // Overall Acceleration of the center of mass in X direction
        customvalue[31] = (customvalue[24] + par[3]*par[10]*(customvalue[26]*cos(x[4]) - x[5]*x[5]*sin(x[4])))/(1 + par[3]);
        // Overall Acceleration of the center of mass in Z direction
        customvalue[32] = (customvalue[25] - par[3]*par[10]*(customvalue[26]*sin(x[4]) + x[5]*x[5]*cos(x[4])))/(1 + par[3]);
        // Overall Input Base Excitation Displacement in X direction
        customvalue[33] = rb*sin((pi/180)*par[2]);
        // Overall Input Base Excitation Displacement in Z direction
        customvalue[34] = rb*cos((pi/180)*par[2]);
        // Overall Input Base Excitation Velocity in X direction
        customvalue[35] = drb*sin((pi/180)*par[2]);
        // Overall Input Base Excitation Velocity in Z direction
        customvalue[36] = drb*cos((pi/180)*par[2]);
        // Overall Input Base Excitation Acceleration in X direction
        customvalue[37] = ddrb*sin((pi/180)*par[2]);
        // Overall Input Base Excitation Acceleration in Z direction
        customvalue[38] = ddrb*cos((pi/180)*par[2]);
        // Overall minimum value of the X acceleration of the system
        min_value(customvalue[24], &customvalue[39]);
        // Overall maximum value of the X acceleration of the system
        max_value(customvalue[24], &customvalue[40]);
        // Overall minimum value of the Z acceleration of the system
        min_value(customvalue[25], &customvalue[41]);
        // Overall maximum value of the Z acceleration of the system
        max_value(customvalue[25], &customvalue[42]);
        // Overall minimum value of the Phi acceleration of the system
        min_value(customvalue[26], &customvalue[43]);
        // Overall maximum value of the Phi acceleration of the system
        min_value(customvalue[26], &customvalue[44]);
        // Overall minimum value of the displacement of the center of mass in X direction
        // Overall maximum value of the displacement of the center of mass in Z direction
        // Overall minimum value of the velocity of the center of mass in X direction
        // Overall minimum value of the velocity of the center of mass in Z direction
        // Overall minimum value of the acceleration of the center of mass in X direction
        // Overall minimum value of the acceleration of the center of mass in Z direction
        // Accumulate the value of the square of the overall X acceleration of the system
        customvalue[45] = RMS(&customvalue[45], customvalue[24], N, 0);
        // Accumulate the valur of the square of the overall Z acceleration of the system
        customvalue[46] = RMS(&customvalue[46], customvalue[25], N, 0);
        // Accumulate the value of the square of the overall Phi acceleration of the system
        customvalue[47] = RMS(&customvalue[47], customvalue[26], N, 0);
        // Accumulate the value of the square of the Overall displacement of the center of mass in X direction
        customvalue[48] = RMS(&customvalue[48], customvalue[27], N, 0);
        // Accumulate the value of the square of the Overall displacement of the center of mass in Z direction
        customvalue[49] = RMS(&customvalue[49], customvalue[28], N, 0);
        // Accumulate the value of the square of the Overall velocity of the center of mass in X direction
        customvalue[50] = RMS(&customvalue[50], customvalue[29], N, 0);
        // Accumulate the value of the square of the Overall velocity of the center of mass in Z direction
        customvalue[51] = RMS(&customvalue[51], customvalue[30], N, 0);
        // Accumulate the value of the square of the Overall acceleration of the center of mass in X direction
        customvalue[52] = RMS(&customvalue[52], customvalue[31], N, 0);
        // Accumulate the value of the square of the Overall acceleration of the center of mass in Z direction
        customvalue[53] = RMS(&customvalue[53], customvalue[32], N, 0);
        // Accumulate the value of the square of the overall input base excitation displacement in X direction
        customvalue[54] = RMS(&customvalue[54], customvalue[33], N, 0);
        // Accumulate the value of the square of the overall input base excitation displacement in Z direction
        customvalue[55] = RMS(&customvalue[55], customvalue[34], N, 0);
        // Accumulate the value of the square of the overall input base excitation velocity in X direction
        customvalue[56] = RMS(&customvalue[56], customvalue[35], N, 0);
        // Accumulate the value of the square of the overall input base excitation velocity in Z direction
        customvalue[57] = RMS(&customvalue[57], customvalue[36], N, 0);
        // Accumulate the value of the square of the overall input base excitation acceleration in X direction
        customvalue[58] = RMS(&customvalue[58], customvalue[37], N, 0);
        // Accumulate the value of the square of the overall input base excitation acceleration in Z direction
        customvalue[59] = RMS(&customvalue[59], customvalue[38], N, 0);
    } 
    // Mode to perform calculations at the end of the time series
    else if (mode == 3) {
        // RMS X acceleration of the system in the steady state regime
        customvalue[47] = RMS(&customvalue[14], customvalue[0], (int)(N*steadystateperc), 1);
        // RMS Z acceleration of the system in the steady state regime
        customvalue[48] = RMS(&customvalue[15], customvalue[1], (int)(N*steadystateperc), 1);
        // RMS Phi acceleration of the system in the steady state regime
        customvalue[49] = RMS(&customvalue[16], customvalue[2], (int)(N*steadystateperc), 1);
        // Overall RMS X acceleration of the system
        customvalue[50] = RMS(&customvalue[38], customvalue[17], N, 1);
        // Overall RMS Z acceleration of the system
        customvalue[51] = RMS(&customvalue[39], customvalue[18], N, 1);
        // Overall RMS Phi acceleration of the system
        customvalue[52] = RMS(&customvalue[40], customvalue[19], N, 1);
        // Overall RMS of the input base excitation displacement in X direction
        customvalue[53] = RMS(&customvalue[41], customvalue[26], N, 1);
        // Overall RMS of the input base excitation displacement in Z direction
        customvalue[54] = RMS(&customvalue[42], customvalue[27], N, 1);
        // Overall RMS of the input base excitation velocity in X direction
        customvalue[55] = RMS(&customvalue[43], customvalue[28], N, 1);
        // Overall RMS of the input base excitation velocity in Z direction
        customvalue[56] = RMS(&customvalue[44], customvalue[29], N, 1);
        // Overall RMS of the input base excitation acceleration in X direction
        customvalue[57] = RMS(&customvalue[45], customvalue[30], N, 1);
        // Overall RMS of the input base excitation acceleration in Z direction
        customvalue[58] = RMS(&customvalue[46], customvalue[31], N, 1);
        // Peak to Peak Average Electrical Output Power of the Piezoelectric Element
        customvalue[59] = (par[15]*par[11]/par[12])*xrms[6];
        // Peak to Peak Average Electrical Output Power of the Electromagnetic Converter
        customvalue[60] = (par[16]*par[13]/par[14])*xrms[7];
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