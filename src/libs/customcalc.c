#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nlosc.h"
#include "nldyn.h"
#include "odesystems.h"
#include "defines.h"
#include "msg.h"

// Methods
static void assign_names(char **strings, const int nvalues, char **names) {
    // Get every row of strings and copy to names
    for (int i = 0; i < nvalues; i++) {
        if (strlen(strings[i]) + 1 >= MAX_CCALC_NAME_LEN) {
            print_error("CUSTOM CALCULATIONS FUNCTION ERROR: One or more given names of the custom values are too big, please assign smaller names of maximum length of %d before trying to run the program.\n", MAX_CCALC_NAME_LEN);
            print_exit_prog();
            exit(EXIT_FAILURE);
        }
        strcpy(names[i], strings[i]);
    }
}

static void error(int mode) {
    print_debug("Custom Function using mode = %d, please use 0, 1, 2 or 3", mode);
    print_exit_prog();
    exit(EXIT_FAILURE);
}

static void spins(double initial_angle, double *previous_angle, double *current_angle, double angle, double *positive_spin, double *negative_spin, int index) {
    // Store previous and current angle
    if (index == 0) {
        (*previous_angle) = initial_angle;
        (*current_angle) = angle;
    }
    else {
        (*previous_angle) = (*current_angle);
        (*current_angle) = angle;
    }
    // Check if the current angle is smaller or bigger than the previous angle and account the advance or retreat of the angle
    if ((*current_angle) > (*previous_angle)) {
        (*positive_spin) = (*positive_spin) + fabs((*current_angle));
    }
    else if ((*current_angle) < (*previous_angle)) {
        (*negative_spin) = (*negative_spin) + fabs((*current_angle));
    }
}

static void time_to_flip(double t, double initial_angle, double current_angle, double *tflip) {
    if ((*tflip) == 0.0) {
        const double pi = 4 * atan(1);  // Pi number definition
        // If flipping happens, store the time and mark as done
        if ((current_angle >= initial_angle + 2*pi) || (current_angle <= initial_angle - 2*pi)) {
            (*tflip) = t;
        }
    }
    else {
        return;
    }
}

// Methods for pend_oscillator_EH
static double pend_oscillator_XCM(double X, double rho, double l, double Phi) {
    return ((X + rho*l*sin(Phi))/(1 + rho)); 
}

static double pend_oscillator_ZCM(double Z, double rho, double l, double Phi) {
    return (Z + rho*l*cos(Phi)/(1 + rho));
}

static double pend_oscillator_dXCM(double dX, double rho, double l, double Phi, double dPhi) {
    return ((dX + rho*l*dPhi*cos(Phi))/(1 + rho));
}

static double pend_oscillator_dZCM(double dZ, double rho, double l, double Phi, double dPhi) {
    return ((dZ - rho*l*dPhi*sin(Phi))/(1 + rho));
}

static double pend_oscillator_ddXCM(double ddX, double rho, double l, double Phi, double dPhi, double ddPhi) {
    return ((ddX + rho*l*(ddPhi*cos(Phi) - dPhi*dPhi*sin(Phi)))/(1 + rho));
}

static double pend_oscillator_ddZCM(double ddZ, double rho, double l, double Phi, double dPhi, double ddPhi) {
    return ((ddZ - rho*l*(ddPhi*sin(Phi) + dPhi*dPhi*cos(Phi)))/(1 + rho));
}


// Methods for duffing_2DoF_EH
static double duffing_2DoF_EH_XCM(double X1, double X2, double rho) {
    return ((X1 + rho*X2)/(1 + rho));
}

static double duffing_2DoF_EH_dXCM(double dX1, double dX2, double rho) {
    return ((dX1 + rho*dX2)/(1 + rho));
}

static double duffing_2DoF_EH_ddXCM(double ddX1, double ddX2, double rho) {
    return ((ddX1 + rho*ddX2)/(1 + rho));
}


// Custom Calculations
void customcalc(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
    return;
}

void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
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
                            "TotalPout"
                        };
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames);
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
        //customvalue[19] = (par[5]*par[6]/par[7])*xrms[2]*xrms[2];
        customvalue[19] = par[6]*xrms[2]*xrms[2];
    }
    else {
        error(mode);
    }
}

void customcalc_pend_oscillator_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
    /* OMEGA   = par[0]   |   zeta_z    = par[5]   |   l         = par[10]   |   chi_PZ = par[15]       |   x[0] = x       |   x[5] = dphi/dt
       gamma   = par[1]   |   zeta_t    = par[6]   |   varphi_PZ = par[11]   |   chi_EM = par[16]       |   x[1] = dx/dt   |   x[6] = v
       mu      = par[2]   |   OMEGA_s   = par[7]   |   kappa_PZ  = par[12]   |                          |   x[2] = z       |   x[7] = i
       rho     = par[3]   |   OMEGA_phi = par[8]   |   varphi_EM = par[13]   |                          |   x[3] = dz/dt   |
       zeta_x  = par[4]   |   OMEGA_t   = par[9]   |   kappa_EM  = par[14]   |                          |   x[4] = phi     |                   */
    
    const double pi = 4 * atan(1);
    double rb = par[1]*sin(par[0]*t);
    double drb = par[1]*par[0]*cos(par[0]*t);
    double ddrb = -par[1]*par[0]*par[0]*sin(par[0]*t);

    double f[8];
    pend_oscillator_EH(8, x, t, par, f);
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
                            "PoutPZ_Avg", "PoutEM_Avg",
                            "prev_ang", "curnt_ang", "pos_spin", "neg_spin",
                            "OVRLL_prev_ang", "OVRLL_curnt_ang", "OVRLL_pos_spin", "OVRLL_neg_spin",
                            "tflip"
                            };
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames);
    }
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) { 
        // X Acceleration of the system in steady state regime ("ddX")
        customvalue[0] = f[1];
        // Z Acceleration of the system in steady state regime ("ddZ")
        customvalue[1] = f[3];
        // Phi Acceleration of the system in steady state regime ("ddPhi")
        customvalue[2] = f[5];
        // Displacement of the center of mass in X direction ("Xcm")
        customvalue[3] = pend_oscillator_XCM(x[0], par[3], par[10], x[4]);          // (x[0] + par[3]*par[10]*sin(x[4]))/(1 + par[3]);
        // Displacement of the center of mass in Z direction ("Zcm")
        customvalue[4] = pend_oscillator_ZCM(x[2], par[3], par[10], x[4]);          // (x[2] + par[3]*par[10]*cos(x[4]))/(1 + par[3]);
        // Velocity of the center of mass in X direction ("dXcm")
        customvalue[5] = pend_oscillator_dXCM(x[1], par[3], par[10], x[4], x[5]);   // (x[1] + par[3]*par[10]*x[5]*cos(x[4]))/(1 + par[3]);
        // Velocity of the center of mass in Z direction ("dZcm")
        customvalue[6] = pend_oscillator_dZCM(x[3], par[3], par[10], x[4], x[5]);   // (x[3] - par[3]*par[10]*x[5]*sin(x[4]))/(1 + par[3]);
        // Acceleration of the center of mass in X direction ("ddXcm")
        customvalue[7] = pend_oscillator_ddXCM(f[1], par[3], par[10], x[4], x[5], f[5]);  //(customvalue[0] + par[3]*par[10]* (customvalue[2]*cos(x[4]) - x[5]*x[5]*sin(x[4])))/(1 + par[3]);
        // Acceleration of the center of mass in Z direction ("ddZcm")
        customvalue[8] = pend_oscillator_ddZCM(f[3], par[3], par[10], x[4], x[5], f[5]);  //(customvalue[1] - par[3]*par[10]* (customvalue[2]*sin(x[4]) + x[5]*x[5]*cos(x[4])))/(1 + par[3]);
        // Minimum values of the accelerations in steady state regime ("ddX_MIN", "ddZ_MIN", "ddPhi_MIN")
        for (int i = 0; i < 3; i++) { // From customvalue[9] to customvalue[11]
            if (currenttimestep == N*steadystateperc) {
                customvalue[9+i] = customvalue[i];    
            }
            else {
                min_value(customvalue[i], &customvalue[9+i]);
            }
        }
        // Maximum values of the accelerations in steady state regime ("ddX_MAX", "ddZ_MAX", "ddPhi_MAX")
        for (int i = 0; i < 3; i++) { // From customvalue[12] to customvalue[14]
            if (currenttimestep == N*steadystateperc) {
                customvalue[12+i] = customvalue[i];
            }
            else {
                max_value(customvalue[i], &customvalue[12+i]);
            }            
        }
        // Minimum values of the values of the center of mass in steady state regime ("Xcm_MIN", "Zcm_MIN", "dXcm_MIN", "dZcm_MIN", "ddXcm_MIN", "ddZcm_MIN")
        for (int i = 0; i < 6; i++) { // From customvalue[15] to customvalue[20]
            if (currenttimestep == N*steadystateperc) {
                customvalue[15+i] = customvalue[3+i];
            }
            else {
                min_value(customvalue[3+i], &customvalue[15+i]);
            }
        }
        // Maximum values of the values of the center of mass in steady state regime ("Xcm_MAX", "Zcm_MAX", "dXcm_MAX", "dZcm_MAX", "ddXcm_MAX", "ddZcm_MAX")
        for (int i = 0; i < 6; i++) { // From customvalue[21] to customvalue[26]
            if (currenttimestep == N*steadystateperc) {
                customvalue[21+i] = customvalue[3+i];
            }
            else {
                max_value(customvalue[3+i], &customvalue[21+i]);
            }
        }
        // Accumulate the value of the square of accelerations of the system in steady state regime ("Sum(ddX^2)", "Sum(ddZ^2)", "Sum(ddPhi^2)")
        for (int i = 0; i < 3; i++) { // From customvalue[27] to customvalue[29]
            customvalue[27+i] = RMS(&customvalue[27+i], customvalue[i], (int)(N*steadystateperc), 0);
        }
        // Accumulate the value of the square of values of the center of mass in steady state regime ("Sum(Xcm^2)", "Sum(Zcm^2)", "Sum(dXcm^2)", "Sum(dZcm^2)", "Sum(ddXcm^2)", "Sum(ddZcm^2)")
        for (int i = 0; i < 6; i++) { // From customvalue[30] to customvalue[35]
            customvalue[30+i] = RMS(&customvalue[30+i], customvalue[3+i], (int)(N*steadystateperc), 0);
        }     

        /* Compute positive and negative spinning of the pendulum in steady state regime */
        // "prev_ang", "curnt_ang", "pos_spin", "neg_spin"
        spins(x[4], &customvalue[110], &customvalue[111], x[4], &customvalue[112], &customvalue[113], currenttimestep);
    }
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        // Overall X Acceleration of the system ("OVRLL_ddX")
        customvalue[36] = f[1];
        // Overall Z Acceleration of the system ("OVRLL_ddZ")
        customvalue[37] = f[3];
        // Overall Phi Acceleration of the system ("OVRLL_ddPhi")
        customvalue[38] = f[5];
        // Overall displacement of the center of mass in X direction ("OVRLL_Xcm")
        customvalue[39] = pend_oscillator_XCM(x[0], par[3], par[10], x[4]);            //(x[0] + par[3]*par[10]*sin(x[4]))/(1 + par[3]); 
        // Overall displacement of the center of mass in Z direction ("OVRLL_Zcm")
        customvalue[40] = pend_oscillator_ZCM(x[2], par[3], par[10], x[4]);            //(x[2] + par[3]*par[10]*cos(x[4]))/(1 + par[3]);
        // Overall velocity of the center of mass in X direction ("OVRLL_dXcm")
        customvalue[41] = pend_oscillator_dXCM(x[1], par[3], par[10], x[4], x[5]);     //(x[1] + par[3]*par[10]*x[5]*cos(x[4]))/(1 + par[3]);
        // Overall velocity of the center of mass in Z direction ("OVRLL_dZcm")
        customvalue[42] = pend_oscillator_dZCM(x[3], par[3], par[10], x[4], x[5]);     //(x[3] - par[3]*par[10]*x[5]*sin(x[4]))/(1 + par[3]);
        // Overall Acceleration of the center of mass in X direction ("OVRLL_ddXcm")
        customvalue[43] = pend_oscillator_ddXCM(customvalue[36], par[3], par[10], x[4], x[5], customvalue[38]); //(customvalue[36] + par[3]*par[10]*(customvalue[38]*cos(x[4]) - x[5]*x[5]*sin(x[4])))/(1 + par[3]);
        // Overall Acceleration of the center of mass in Z direction ("OVRLL_ddZcm")
        customvalue[44] = pend_oscillator_ddZCM(customvalue[37], par[3], par[10], x[4], x[5], customvalue[38]);   //(customvalue[37] - par[3]*par[10]*(customvalue[38]*sin(x[4]) + x[5]*x[5]*cos(x[4])))/(1 + par[3]);
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
        /* Overall Minimum values of the values of the center of mass ("OVRLL_Xcm_MIN", "OVRLL_Zcm_MIN", "OVRLL_dXcm_MIN", "OVRLL_dZcm_MIN", "OVRLL_ddXcm_MIN", "OVRLL_ddZcm_MIN") */
        if (currenttimestep == 0) {  // Identify initial conditions of the center of mass variables
            customvalue[57] = pend_oscillator_XCM(IC[0], par[3], par[10], IC[4]);
            customvalue[58] = pend_oscillator_ZCM(IC[2], par[3], par[10], IC[4]);
            customvalue[59] = pend_oscillator_dXCM(IC[1], par[3], par[10], IC[4], IC[5]);
            customvalue[60] = pend_oscillator_dZCM(IC[3], par[3], par[10], IC[4], IC[5]);
            // customvalue[61] and customvalue[62] are already zero and don't require calculations of initial condition
        }
        for (int i = 0; i < 6; i++) { // From customvalue[57] to customvalue[62]
            min_value(customvalue[39+i], &customvalue[57+i]);
        }
        /* Overall Maximum values of the values of the center of mass ("OVRLL_Xcm_MAX", "OVRLL_Zcm_MAX", "OVRLL_dXcm_MAX", "OVRLL_dZcm_MAX", "OVRLL_ddXcm_MAX", "OVRLL_ddZcm_MAX") */
        if (currenttimestep == 0) {  // Identify initial conditions of the center of mass variables
            customvalue[63] = pend_oscillator_XCM(IC[0], par[3], par[10], IC[4]);
            customvalue[64] = pend_oscillator_ZCM(IC[2], par[3], par[10], IC[4]);
            customvalue[65] = pend_oscillator_dXCM(IC[1], par[3], par[10], IC[4], IC[5]);
            customvalue[66] = pend_oscillator_dZCM(IC[3], par[3], par[10], IC[4], IC[5]);
            // customvalue[67] and customvalue[68] are already zero and don't require calculations of initial conditions
        }
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
        
        /* Compute positive and negative spinning of the pendulum */
        // "OVRLL_prev_ang", "OVRLL_curnt_ang", "OVRLL_pos_spin", "OVRLL_neg_spin"
        spins(IC[4], &customvalue[114], &customvalue[115], x[4], &customvalue[116], &customvalue[117], currenttimestep);

        /* Time To Flip Calculations */ 
        // "tflip"
        time_to_flip(t, IC[4], x[4], &customvalue[118]);    
    } 
    // Mode to perform calculations at the end of the time series    
    else if (mode == 3) {
        // RMS acceleration of the system in steady state regime ("ddX_RMS", "ddZ_RMS", "ddPhi_RMS")
        for (int i = 0; i < 3; i++) { // From customvalue[84] to customvalue[86]
            customvalue[84+i] = RMS(&customvalue[27+i], customvalue[i], (int)(N*steadystateperc), 1);
        }
        // RMS of the center of mass in steady state regime ("Xcm_RMS", "Zcm_RMS", "dXcm_RMS", "dZcm_RMS", "ddXcm_RMS", "ddZcm_RMS")
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
        //customvalue[108] = (par[15]*par[11]/par[12])*xrms[6]*xrms[6];
        customvalue[108] = par[11]*xrms[6]*xrms[6];
        // Peak to Peak Average Electrical Output Power of the Electromagnetic Converter
        //customvalue[109] = (par[16]*par[13]/par[14])*xrms[7]*xrms[7];
        customvalue[109] = par[13]*xrms[7]*xrms[7];
    }
    else {
        error(mode);
    }
}

void customcalc_duffing_2DoF_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
    /* OMEGA   = par[0]   |   alpha_1 = par[5]   |   chi_1    = par[10]   |   kappa_2 = par[15]   |   x1  = x[0] | v2 = x[5]
       gamma   = par[1]   |   alpha_2 = par[6]   |   chi_2    = par[11]   |                       |   dx1 = x[1] |    
       rho     = par[2]   |   beta_1  = par[7]   |   varphi_1 = par[12]   |                       |   x2  = x[2] |
       zeta_1  = par[3]   |   beta_2  = par[8]   |   varphi_2 = par[13]   |                       |   dx2 = x[3] |
       zeta_2  = par[4]   |   OMEGA_s = par[9]   |   kappa_1  = par[14]   |                       |   v1  = x[4] |           */
    
    double f[6];
    duffing_2DoF_EH(6, x, t, par, f);
    
    // Mode to define names to be printed in the output file
    if (mode == 0) {    
        char *names[] = {   "ddX1", "ddX2", "Xcm", "dXcm", "ddXcm", 
                            "ddX1_MIN", "ddX2_MIN",
                            "ddX1_MAX", "ddX2_MAX",
                            "Xcm_MIN", "dXcm_MIN", "ddXcm_MIN",
                            "Xcm_MAX", "dXcm_MAX", "ddXcm_MAX",
                            "Sum(ddX1^2)", "Sum(ddX2^2)", 
                            "Sum(Xcm^2)", "Sum(dXcm^2)", "Sum(ddXcm^2)",    
                            "OVRL_ddX1", "OVRL_ddX2", "OVRL_Xcm", "OVRL_dXcm", "OVRL_ddXcm", "OVRL_Xb", "OVRL_dXb", "OVRL_ddXb",
                            "OVRL_ddX1_MIN", "OVRL_ddX2_MIN", 
                            "OVRL_ddX1_MAX", "OVRL_ddX2_MAX",
                            "OVRL_Xcm_MIN", "OVRL_dXcm_MIN", "OVRL_ddXcm_MIN",
                            "OVRL_Xcm_MAX", "OVRL_dXcm_MAX", "OVRL_ddXcm_MAX",
                            "Sum(OVRL_ddX1^2)", "Sum(OVRL_ddX2^2)", 
                            "Sum(OVRL_Xcm^2)", "Sum(OVRL_dXcm^2)", "Sum(OVRL_ddXcm^2)",
                            "Sum(OVRL_Xb)", "Sum(OVRL_dXb)", "Sum(OVRL_ddXb)",
                            "ddX1_RMS", "ddX2_RMS",
                            "Xcm_RMS", "dXcm_RMS", "ddXcm_RMS",
                            "OVRL_ddX1_RMS", "OVRL_ddX2_RMS",
                            "OVRL_Xcm_RMS", "OVRL_dXcm_RMS", "OVRL_ddXcm_RMS",
                            "OVRL_Xb_RMS", "OVRL_dXb_RMS", "OVRL_ddXb_RMS",
                            "Pout1", "Pout2", "TotalPout", "PoutDens",
                            "Xrel", "dXrel", "ddXrel", 
                            "Xrel_MIN", "dXrel_MIN", "ddXrel_MIN",
                            "Xrel_MAX", "dXrel_MAX", "ddXrel_MAX",
                            "Sum(Xrel^2)", "Sum(dXrel^2)", "Sum(ddXrel^2)",
                            "OVRL_Xrel", "OVRL_dXrel", "OVRL_ddXrel", 
                            "OVRL_Xrel_MIN", "OVRL_dXrel_MIN", "OVRL_ddXrel_MIN",
                            "OVRL_Xrel_MAX", "OVRL_dXrel_MAX", "OVRL_ddXrel_MAX",
                            "Sum(OVRL_Xrel^2)", "Sum(OVRL_dXrel^2)", "Sum(OVRL_ddXrel^2)",
                            "Xrel_RMS", "dXrel_RMS", "ddXrel_RMS",
                            "OVRL_Xrel_RMS", "OVRL_dXrel_RMS", "OVRL_ddXrel_RMS", 
                        };
        
        
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames);
    } 
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) {
        // X1 Acceleration of the system in steady state regime ("ddX1")
        customvalue[0] = f[1];
        // X2 Acceleration of the system in steady state regime ("ddX2")
        customvalue[1] = f[3];
        // Displacement of the center of mass ("Xcm")
        customvalue[2] = duffing_2DoF_EH_XCM(x[0], x[2], par[2]);
        // Velocity of the center of mass ("dXcm")
        customvalue[3] = duffing_2DoF_EH_dXCM(x[1], x[3], par[2]);
        // Acceleration of the center of mass ("ddXcm")
        customvalue[4] = duffing_2DoF_EH_ddXCM(f[1], f[3], par[2]);
        
        // X2 - X1 Displacement of the system in steady state regime ("Xrel")
        customvalue[63] = x[2] - x[0];
        // dX2 - dX1 Velocity of the system in steady state regime ("dXrel")
        customvalue[64] = x[3] - x[1];
        // ddX2 - ddX1 Acceleration of the system in steady state regime ("ddXrel")
        customvalue[65] = f[3] - f[1];

        for (int i = 0; i < 2; i++) { 
            // Minimum and Maximum values of the accelerations in steady state regime ("ddX1_MIN", "ddX2_MIN", "ddX1_MAX", "ddX2_MAX")
            if (currenttimestep == N*steadystateperc) {
                customvalue[5+i] = customvalue[i]; // From customvalue[5] to customvalue[6]    
                customvalue[7+i] = customvalue[i]; // From customvalue[7] to customvalue[8]
            }
            else {
                min_value(customvalue[i], &customvalue[5+i]); // "ddX1_MIN", "ddX2_MIN"
                max_value(customvalue[i], &customvalue[7+i]); // "ddX1_MAX", "ddX2_MAX"
            }
            // Accumulate the value of the square of accelerations of the system in steady state regime ("Sum(ddX1^2)", "Sum(ddX2^2)")
            customvalue[15+i] = RMS(&customvalue[15+i], customvalue[i], (int)(N*steadystateperc), 0);  // From customvalue[15] to customvalue[16] 
        }
           
        for (int i = 0; i < 3; i++) {   
            // Minimum and Maximum values of the values in steady state regime  
            if (currenttimestep == N*steadystateperc) {
                customvalue[9+i] = customvalue[2+i];   // From customvalue[9] to customvalue[11] 
                customvalue[12+i] = customvalue[2+i];  // From customvalue[12] to customvalue[14] 
                customvalue[66+i] = customvalue[63+i]; // From customvalue[66] to customvalue[68] 
                customvalue[69+i] = customvalue[63+i]; // From customvalue[69] to customvalue[71] 
            }
            else {
                min_value(customvalue[2+i], &customvalue[9+i]);  // "Xcm_MIN", "dXcm_MIN", "ddXcm_MIN"
                max_value(customvalue[2+i], &customvalue[12+i]); // "Xcm_MAX", "dXcm_MAX", "ddXcm_MAX"
                min_value(customvalue[63+i], &customvalue[66+i]); // "Xrel_MIN", "dXrel_MIN", "ddXrel_MIN"
                max_value(customvalue[63+i], &customvalue[69+i]); // "Xrel_MAX", "dXrel_MAX", "ddXrel_MAX"
            }               
            // Accumulate the value of the square of values of the center of mass in steady state regime ("Sum(Xcm^2)", "Sum(dXcm^2)", "Sum(ddXcm^2)")
            customvalue[17+i] = RMS(&customvalue[17+i], customvalue[2+i], (int)(N*steadystateperc), 0); // From customvalue[17] to customvalue[19] 
            // Accumulate the value of the square of values of the relative motion X2-X1 and its derivatives in steady state regime ("Sum(Xrel^2)", "Sum(dXrel^2)", "Sum(ddXrel^2)")
            customvalue[72+i] = RMS(&customvalue[72+i], customvalue[63+i], (int)(N*steadystateperc), 0); // From customvalue[72] to customvalue[74]    
        }
    } 
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        // Overall X Acceleration of the system ("OVRLL_ddX1")
        customvalue[20] = f[1];
        // Overall Z Acceleration of the system ("OVRLL_ddX2")
        customvalue[21] = f[3];
        // Overall Displacement of the center of mass ("OVRLL_Xcm")
        customvalue[22] = duffing_2DoF_EH_XCM(x[0], x[2], par[2]);
        // Overall Velocity of the center of mass ("OVRLL_dXcm")
        customvalue[23] = duffing_2DoF_EH_dXCM(x[1], x[3], par[2]);
        // Overall Acceleration of the center of mass ("OVRLL_ddXcm")
        customvalue[24] = duffing_2DoF_EH_ddXCM(f[1], f[3], par[2]);
        // Overall Input Base Excitation Displacement ("OVRLL_Xb")    
        customvalue[25] = par[1]*sin(par[0]*t);
        // Overall Input Base Excitation Velocity ("OVRLL_dXb")
        customvalue[26] = par[1]*par[0]*cos(par[0]*t);
        // Overall Input Base Excitation Acceleration ("OVRLL_ddXb")
        customvalue[27] = -par[1]*par[0]*par[0]*sin(par[0]*t);
        
        // Overall X2 - X1 Displacement of the system in steady state regime ("OVRL_Xrel")
        customvalue[75] = x[2] - x[0];
        // Overall dX2 - dX1 Velocity of the system in steady state regime ("OVRL_dXrel")
        customvalue[76] = x[3] - x[1];
        // Overall ddX2 - ddX1 Acceleration of the system in steady state regime ("OVRL_ddXrel")
        customvalue[77] = f[3] - f[1];

        for (int i = 0; i < 2; i++) { 
            // Overall Minimum and Maximum values of the accelerations ("OVRL_ddX1_MIN", "OVRL_ddX2_MIN", "OVRL_ddX1_MAX", "OVRL_ddX2_MAX")
            min_value(customvalue[20+i], &customvalue[28+i]);   // From customvalue[28] to customvalue[29]
            max_value(customvalue[20+i], &customvalue[30+i]);   // From customvalue[30] to customvalue[31]
            // Accumulate the overall value of the square of accelerations of the system ("OVRL_Sum(ddX1^2)", "OVRL_Sum(ddX2^2)")
            customvalue[38+i] = RMS(&customvalue[38+i], customvalue[20+i], N, 0); // From customvalue[38] to customvalue[39]
        }

        for (int i = 0; i < 3; i++) {
            // Overall Minimum and Maximum values 
            min_value(customvalue[22+i], &customvalue[32+i]);   // From customvalue[32] to customvalue[34] ("OVRL_Xcm_MIN", "OVRL_dXcm_MIN", "OVRL_ddXcm_MIN")
            max_value(customvalue[22+i], &customvalue[35+i]);   // From customvalue[35] to customvalue[37] ("OVRL_Xcm_MAX", "OVRL_dXcm_MAX", "OVRL_ddXcm_MAX")
            min_value(customvalue[75+i], &customvalue[78+i]);   // From customvalue[78] to customvalue[80] ("OVRL_Xrel_MIN", "OVRL_dXrel_MIN", "OVRL_ddXrel_MIN")
            max_value(customvalue[75+i], &customvalue[81+i]);   // From customvalue[81] to customvalue[83] ("OVRL_Xrel_MAX", "OVRL_dXrel_MAX", "OVRL_ddXrel_MAX")
            // Accumulate the overall value of the square of values of the center of mass ("OVRL_Sum(Xcm^2)", "OVRL_Sum(dXcm^2)", "OVRL_Sum(ddXcm^2)")
            customvalue[40+i] = RMS(&customvalue[40+i], customvalue[22+i], N, 0); // From customvalue[40] to customvalue[42]
            // Accumulate the value of the square of the overall input base excitation values ("OVRL_Sum(Xb^2)", "OVRL_Sum(dXb^2)", "OVRL_Sum(ddXb^2)")
            customvalue[43+i] = RMS(&customvalue[43+i], customvalue[25+i], N, 0); // From customvalue[43] to customvalue[45]
            // Accumulate the value of the square of the overall Xrel=X2-X1 values and its derivatives ("OVRL_Sum(Xrel^2)", "OVRL_Sum(dXrel^2)", "OVRL_Sum(ddXrel^2)")
            customvalue[84+i] = RMS(&customvalue[84+i], customvalue[75+i], N, 0); // From customvalue[84] to customvalue[86]
        }
    }
    // Mode to perform calculations at the end of the time series    
    else if (mode == 3) {

        for (int i = 0; i < 2; i++) { 
            // RMS acceleration of the system in steady state regime ("ddX1_RMS", "ddX2_RMS")
            customvalue[46+i] = RMS(&customvalue[15+i], customvalue[i], (int)(N*steadystateperc), 1); // From customvalue[46] to customvalue[47]
            // Overall RMS acceleration of the system in steady state regime ("OVRLL_ddX1_RMS", "OVRLL_ddX2_RMS")
            customvalue[51+i] = RMS(&customvalue[38+i], customvalue[20+i], N, 1); // From customvalue[51] to customvalue[52]
        }
        
        for (int i = 0; i < 3; i ++) {
            // RMS of the center of mass in steady state regime ("Xcm_RMS", "dXcm_RMS", "ddXcm_RMS")
            customvalue[48+i] = RMS(&customvalue[17+i], customvalue[2+i], (int)(N*steadystateperc), 1); // From customvalue[48] to customvalue[50]
            // Overall RMS of the center of mass in steady state regime ("OVRLL_Xcm_RMS", "OVRLL_dXcm_RMS", "OVRLL_ddXcm_RMS")
            customvalue[53+i] = RMS(&customvalue[40+i], customvalue[22+i], N, 1); // From customvalue[53] to customvalue[55]
            // Overal RMS of the input base excitation values ("OVRLL_Xb_RMS", "OVRLL_dXb_RMS", "OVRLL_ddXb_RMS",)
            customvalue[56+i] = RMS(&customvalue[43+i], customvalue[25+i], N, 1); // From customvalue[56] to customvalue[58]
            // RMS of the motion of Xrel=X2-X1 and its derivatives ("Xrel_RMS", "dXrel_RMS", "ddXrel_RMS")
            customvalue[87+i] = RMS(&customvalue[72+i], customvalue[63+i], (int)(N*steadystateperc), 1); // From customvalue[87] to customvalue[89]
            // RMS of the motion of Xrel=X2-X1 and its derivatives ("OVRL_Xrel_RMS", "OVRL_dXrel_RMS", "OVRL_ddXrel_RMS")
            customvalue[90+i] = RMS(&customvalue[84+i], customvalue[75+i], N, 1); // From customvalue[90] to customvalue[92]
        }

        // Peak to Peak Average Electrical Output Power of the Piezoelectric Element 1
        //customvalue[59] = ((par[10]*par[12])/par[14])*xrms[4]*xrms[4];
        customvalue[59] = par[12]*xrms[4]*xrms[4];
        // Peak to Peak Average Electrical Output Power of the Piezoelectric Element 2
        //customvalue[60] = ((par[11]*par[13])/par[15])*xrms[5]*xrms[5];
        customvalue[60] = par[13]*xrms[5]*xrms[5];
        // Sum of Peak to Peak Average Electrical Output Power of the Piezoelectric Elements
        customvalue[61] = customvalue[59] + customvalue[60];
        // Output Power Density of the system
        customvalue[62] = customvalue[61]/2;

    }
    else {
        error(mode);
    }
}

void customcalc_linear_2DoF_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
    /* OMEGA   = par[0]   |   OMEGA_s  = par[5]   |   kappa_1 = par[10]   |   x1  = x[0] |  v2 = x[5]
       gamma   = par[1]   |   chi_1    = par[6]   |   kappa_2 = par[11]   |   dx1 = x[1] |    
       rho     = par[2]   |   chi_2    = par[7]   |                       |   x2  = x[2] |
       zeta_1  = par[3]   |   varphi_1 = par[8]   |                       |   dx2 = x[3] |
       zeta_2  = par[4]   |   varphi_2 = par[9]   |                       |   v1  = x[4] |         */
    
    double f[6];
    linear_2DoF_EH(6, x, t, par, f);
    
    // Mode to define names to be printed in the output file
    if (mode == 0) {    
        char *names[] = {   "ddX1", "ddX2", "Xcm", "dXcm", "ddXcm", 
                            "ddX1_MIN", "ddX2_MIN",
                            "ddX1_MAX", "ddX2_MAX",
                            "Xcm_MIN", "dXcm_MIN", "ddXcm_MIN",
                            "Xcm_MAX", "dXcm_MAX", "ddXcm_MAX",
                            "Sum(ddX1^2)", "Sum(ddX2^2)", 
                            "Sum(Xcm^2)", "Sum(dXcm^2)", "Sum(ddXcm^2)",    
                            "OVRL_ddX1", "OVRL_ddX2", "OVRL_Xcm", "OVRL_dXcm", "OVRL_ddXcm", "OVRL_Xb", "OVRL_dXb", "OVRL_ddXb",
                            "OVRL_ddX1_MIN", "OVRL_ddX2_MIN", 
                            "OVRL_ddX1_MAX", "OVRL_ddX2_MAX",
                            "OVRL_Xcm_MIN", "OVRL_dXcm_MIN", "OVRL_ddXcm_MIN",
                            "OVRL_Xcm_MAX", "OVRL_dXcm_MAX", "OVRL_ddXcm_MAX",
                            "Sum(OVRL_ddX1^2)", "Sum(OVRL_ddX2^2)", 
                            "Sum(OVRL_Xcm^2)", "Sum(OVRL_dXcm^2)", "Sum(OVRL_ddXcm^2)",
                            "Sum(OVRL_Xb)", "Sum(OVRL_dXb)", "Sum(OVRL_ddXb)",
                            "ddX1_RMS", "ddX2_RMS",
                            "Xcm_RMS", "dXcm_RMS", "ddXcm_RMS",
                            "OVRL_ddX1_RMS", "OVRL_ddX2_RMS",
                            "OVRL_Xcm_RMS", "OVRL_dXcm_RMS", "OVRL_ddXcm_RMS",
                            "OVRL_Xb_RMS", "OVRL_dXb_RMS", "OVRL_ddXb_RMS",
                            "Pout1", "Pout2", "TotalPout", "PoutDens",
                            "Xrel", "dXrel", "ddXrel", 
                            "Xrel_MIN", "dXrel_MIN", "ddXrel_MIN",
                            "Xrel_MAX", "dXrel_MAX", "ddXrel_MAX",
                            "Sum(Xrel^2)", "Sum(dXrel^2)", "Sum(ddXrel^2)",
                            "OVRL_Xrel", "OVRL_dXrel", "OVRL_ddXrel", 
                            "OVRL_Xrel_MIN", "OVRL_dXrel_MIN", "OVRL_ddXrel_MIN",
                            "OVRL_Xrel_MAX", "OVRL_dXrel_MAX", "OVRL_ddXrel_MAX",
                            "Sum(OVRL_Xrel^2)", "Sum(OVRL_dXrel^2)", "Sum(OVRL_ddXrel^2)",
                            "Xrel_RMS", "dXrel_RMS", "ddXrel_RMS",
                            "OVRL_Xrel_RMS", "OVRL_dXrel_RMS", "OVRL_ddXrel_RMS", 
                        };
        
        
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames);
    } 
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) {
        // X1 Acceleration of the system in steady state regime ("ddX1")
        customvalue[0] = f[1];
        // X2 Acceleration of the system in steady state regime ("ddX2")
        customvalue[1] = f[3];
        // Displacement of the center of mass ("Xcm")
        customvalue[2] = duffing_2DoF_EH_XCM(x[0], x[2], par[2]);
        // Velocity of the center of mass ("dXcm")
        customvalue[3] = duffing_2DoF_EH_dXCM(x[1], x[3], par[2]);
        // Acceleration of the center of mass ("ddXcm")
        customvalue[4] = duffing_2DoF_EH_ddXCM(f[1], f[3], par[2]);
        
        // X2 - X1 Displacement of the system in steady state regime ("Xrel")
        customvalue[63] = x[2] - x[0];
        // dX2 - dX1 Velocity of the system in steady state regime ("dXrel")
        customvalue[64] = x[3] - x[1];
        // ddX2 - ddX1 Acceleration of the system in steady state regime ("ddXrel")
        customvalue[65] = f[3] - f[1];

        for (int i = 0; i < 2; i++) { 
            // Minimum and Maximum values of the accelerations in steady state regime ("ddX1_MIN", "ddX2_MIN", "ddX1_MAX", "ddX2_MAX")
            if (currenttimestep == N*steadystateperc) {
                customvalue[5+i] = customvalue[i]; // From customvalue[5] to customvalue[6]    
                customvalue[7+i] = customvalue[i]; // From customvalue[7] to customvalue[8]
            }
            else {
                min_value(customvalue[i], &customvalue[5+i]); // "ddX1_MIN", "ddX2_MIN"
                max_value(customvalue[i], &customvalue[7+i]); // "ddX1_MAX", "ddX2_MAX"
            }
            // Accumulate the value of the square of accelerations of the system in steady state regime ("Sum(ddX1^2)", "Sum(ddX2^2)")
            customvalue[15+i] = RMS(&customvalue[15+i], customvalue[i], (int)(N*steadystateperc), 0);  // From customvalue[15] to customvalue[16] 
        }
           
        for (int i = 0; i < 3; i++) {   
            // Minimum and Maximum values of the values in steady state regime  
            if (currenttimestep == N*steadystateperc) {
                customvalue[9+i] = customvalue[2+i];   // From customvalue[9] to customvalue[11] 
                customvalue[12+i] = customvalue[2+i];  // From customvalue[12] to customvalue[14] 
                customvalue[66+i] = customvalue[63+i]; // From customvalue[66] to customvalue[68] 
                customvalue[69+i] = customvalue[63+i]; // From customvalue[69] to customvalue[71] 
            }
            else {
                min_value(customvalue[2+i], &customvalue[9+i]);  // "Xcm_MIN", "dXcm_MIN", "ddXcm_MIN"
                max_value(customvalue[2+i], &customvalue[12+i]); // "Xcm_MAX", "dXcm_MAX", "ddXcm_MAX"
                min_value(customvalue[63+i], &customvalue[66+i]); // "Xrel_MIN", "dXrel_MIN", "ddXrel_MIN"
                max_value(customvalue[63+i], &customvalue[69+i]); // "Xrel_MAX", "dXrel_MAX", "ddXrel_MAX"
            }               
            // Accumulate the value of the square of values of the center of mass in steady state regime ("Sum(Xcm^2)", "Sum(dXcm^2)", "Sum(ddXcm^2)")
            customvalue[17+i] = RMS(&customvalue[17+i], customvalue[2+i], (int)(N*steadystateperc), 0); // From customvalue[17] to customvalue[19] 
            // Accumulate the value of the square of values of the relative motion X2-X1 and its derivatives in steady state regime ("Sum(Xrel^2)", "Sum(dXrel^2)", "Sum(ddXrel^2)")
            customvalue[72+i] = RMS(&customvalue[72+i], customvalue[63+i], (int)(N*steadystateperc), 0); // From customvalue[72] to customvalue[74]    
        }
    } 
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        // Overall X Acceleration of the system ("OVRLL_ddX1")
        customvalue[20] = f[1];
        // Overall Z Acceleration of the system ("OVRLL_ddX2")
        customvalue[21] = f[3];
        // Overall Displacement of the center of mass ("OVRLL_Xcm")
        customvalue[22] = duffing_2DoF_EH_XCM(x[0], x[2], par[2]);
        // Overall Velocity of the center of mass ("OVRLL_dXcm")
        customvalue[23] = duffing_2DoF_EH_dXCM(x[1], x[3], par[2]);
        // Overall Acceleration of the center of mass ("OVRLL_ddXcm")
        customvalue[24] = duffing_2DoF_EH_ddXCM(f[1], f[3], par[2]);
        // Overall Input Base Excitation Displacement ("OVRLL_Xb")    
        customvalue[25] = par[1]*sin(par[0]*t);
        // Overall Input Base Excitation Velocity ("OVRLL_dXb")
        customvalue[26] = par[1]*par[0]*cos(par[0]*t);
        // Overall Input Base Excitation Acceleration ("OVRLL_ddXb")
        customvalue[27] = -par[1]*par[0]*par[0]*sin(par[0]*t);
        
        // Overall X2 - X1 Displacement of the system in steady state regime ("OVRL_Xrel")
        customvalue[75] = x[2] - x[0];
        // Overall dX2 - dX1 Velocity of the system in steady state regime ("OVRL_dXrel")
        customvalue[76] = x[3] - x[1];
        // Overall ddX2 - ddX1 Acceleration of the system in steady state regime ("OVRL_ddXrel")
        customvalue[77] = f[3] - f[1];

        for (int i = 0; i < 2; i++) { 
            // Overall Minimum and Maximum values of the accelerations ("OVRL_ddX1_MIN", "OVRL_ddX2_MIN", "OVRL_ddX1_MAX", "OVRL_ddX2_MAX")
            min_value(customvalue[20+i], &customvalue[28+i]);   // From customvalue[28] to customvalue[29]
            max_value(customvalue[20+i], &customvalue[30+i]);   // From customvalue[30] to customvalue[31]
            // Accumulate the overall value of the square of accelerations of the system ("OVRL_Sum(ddX1^2)", "OVRL_Sum(ddX2^2)")
            customvalue[38+i] = RMS(&customvalue[38+i], customvalue[20+i], N, 0); // From customvalue[38] to customvalue[39]
        }

        for (int i = 0; i < 3; i++) {
            // Overall Minimum and Maximum values 
            min_value(customvalue[22+i], &customvalue[32+i]);   // From customvalue[32] to customvalue[34] ("OVRL_Xcm_MIN", "OVRL_dXcm_MIN", "OVRL_ddXcm_MIN")
            max_value(customvalue[22+i], &customvalue[35+i]);   // From customvalue[35] to customvalue[37] ("OVRL_Xcm_MAX", "OVRL_dXcm_MAX", "OVRL_ddXcm_MAX")
            min_value(customvalue[75+i], &customvalue[78+i]);   // From customvalue[78] to customvalue[80] ("OVRL_Xrel_MIN", "OVRL_dXrel_MIN", "OVRL_ddXrel_MIN")
            max_value(customvalue[75+i], &customvalue[81+i]);   // From customvalue[81] to customvalue[83] ("OVRL_Xrel_MAX", "OVRL_dXrel_MAX", "OVRL_ddXrel_MAX")
            // Accumulate the overall value of the square of values of the center of mass ("OVRL_Sum(Xcm^2)", "OVRL_Sum(dXcm^2)", "OVRL_Sum(ddXcm^2)")
            customvalue[40+i] = RMS(&customvalue[40+i], customvalue[22+i], N, 0); // From customvalue[40] to customvalue[42]
            // Accumulate the value of the square of the overall input base excitation values ("OVRL_Sum(Xb^2)", "OVRL_Sum(dXb^2)", "OVRL_Sum(ddXb^2)")
            customvalue[43+i] = RMS(&customvalue[43+i], customvalue[25+i], N, 0); // From customvalue[43] to customvalue[45]
            // Accumulate the value of the square of the overall Xrel=X2-X1 values and its derivatives ("OVRL_Sum(Xrel^2)", "OVRL_Sum(dXrel^2)", "OVRL_Sum(ddXrel^2)")
            customvalue[84+i] = RMS(&customvalue[84+i], customvalue[75+i], N, 0); // From customvalue[84] to customvalue[86]
        }
    }
    // Mode to perform calculations at the end of the time series    
    else if (mode == 3) {

        for (int i = 0; i < 2; i++) { 
            // RMS acceleration of the system in steady state regime ("ddX1_RMS", "ddX2_RMS")
            customvalue[46+i] = RMS(&customvalue[15+i], customvalue[i], (int)(N*steadystateperc), 1); // From customvalue[46] to customvalue[47]
            // Overall RMS acceleration of the system in steady state regime ("OVRLL_ddX1_RMS", "OVRLL_ddX2_RMS")
            customvalue[51+i] = RMS(&customvalue[38+i], customvalue[20+i], N, 1); // From customvalue[51] to customvalue[52]
        }
        
        for (int i = 0; i < 3; i ++) {
            // RMS of the center of mass in steady state regime ("Xcm_RMS", "dXcm_RMS", "ddXcm_RMS")
            customvalue[48+i] = RMS(&customvalue[17+i], customvalue[2+i], (int)(N*steadystateperc), 1); // From customvalue[48] to customvalue[50]
            // Overall RMS of the center of mass in steady state regime ("OVRLL_Xcm_RMS", "OVRLL_dXcm_RMS", "OVRLL_ddXcm_RMS")
            customvalue[53+i] = RMS(&customvalue[40+i], customvalue[22+i], N, 1); // From customvalue[53] to customvalue[55]
            // Overal RMS of the input base excitation values ("OVRLL_Xb_RMS", "OVRLL_dXb_RMS", "OVRLL_ddXb_RMS",)
            customvalue[56+i] = RMS(&customvalue[43+i], customvalue[25+i], N, 1); // From customvalue[56] to customvalue[58]
            // RMS of the motion of Xrel=X2-X1 and its derivatives ("Xrel_RMS", "dXrel_RMS", "ddXrel_RMS")
            customvalue[87+i] = RMS(&customvalue[72+i], customvalue[63+i], (int)(N*steadystateperc), 1); // From customvalue[87] to customvalue[89]
            // RMS of the motion of Xrel=X2-X1 and its derivatives ("OVRL_Xrel_RMS", "OVRL_dXrel_RMS", "OVRL_ddXrel_RMS")
            customvalue[90+i] = RMS(&customvalue[84+i], customvalue[75+i], N, 1); // From customvalue[90] to customvalue[92]
        }

        // Peak to Peak Average Electrical Output Power of the Piezoelectric Element 1
        //customvalue[59] = ((par[6]*par[8])/par[10])*xrms[4]*xrms[4];
        customvalue[59] = par[8]*xrms[4]*xrms[4]; 
        // Peak to Peak Average Electrical Output Power of the Piezoelectric Element 2
        //customvalue[60] = ((par[7]*par[9])/par[11])*xrms[5]*xrms[5];
        customvalue[60] = par[9]*xrms[5]*xrms[5];
        // Sum of Peak to Peak Average Electrical Output Power of the Piezoelectric Elements
        customvalue[61] = customvalue[59] + customvalue[60];
        // Output Power Density of the system
        customvalue[62] = customvalue[61]/2;
        
    }
    else {
        error(mode);
    }
}

void customcalc_adeodato_sma_oscillator(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
    /* System Parameters  |                                                              Shape Memory Properties                                                               |   
       -----------------  -   ----------------------------------------------------------------------------------------------------------------------------------------------   -      
       OMEGA   = par[0]   |   sigma_0    = par[6]   |   beta_0    = par[12]   |   alpha  = par[18]  |  s_2      = par[24]   |   load_ant = par[30]   |   T         = par[36]   |   
       gamma   = par[1]   |   sigma_ant  = par[7]   |   beta_ant  = par[13]   |   Ms     = par[19]  |  Ca       = par[25]   |   load     = par[31]   |   t_ant     = par[37]   |      
       c       = par[2]   |   sigma      = par[8]   |   beta      = par[14]   |   Mf     = par[20]  |  Cm       = par[26]   |   Area     = par[32]   |   E_ant     = par[38]   |           
       k       = par[3]   |   strain_0   = par[9]   |   E_0       = par[15]   |   As     = par[21]  |  strain_r = par[27]   |   L_0      = par[33]   |   alpha_ant = par[39]   |
       m       = par[4]   |   strain_ant = par[10]  |   E         = par[16]   |   Af     = par[22]  |  dir_ant  = par[28]   |   Ea       = par[34]   |                         |
       t_0     = par[5]   |   strain     = par[11]  |   alpha_0   = par[17]   |   s_1    = par[23]  |  dir      = par[29]   |   Em       = par[35]   |                         |   */
    // Mode to define names to be printed in the output file
    if (mode == 0) {    
        char *names[] = { "Sigma_0", "Sigma_ant", "Sigma", "Strain_0", "Strain_ant", "Strain",
                          "E_0", "E_ant", "E", "Alpha_0", "Alpha_ant", "Alpha", "Beta_0", "Beta_ant", "Beta", "dir_ant", "dir", "load_ant", "load"};
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames);

        // Assign SMA Initial Conditions Properties
        if (t == t0) {
            // Sigma
            customvalue[0] = par[6];  customvalue[1] = par[7];   customvalue[2] = par[8];
            // Strain
            customvalue[3] = par[9];  customvalue[4] = par[10];  customvalue[5] = par[11];
            // Beta
            customvalue[6] = par[12]; customvalue[7] = par[13];  customvalue[8] = par[14];
            // E
            customvalue[9] = par[15]; customvalue[10] = par[38];  customvalue[11] = par[16];
            // Alpha
            customvalue[12] = par[17]; customvalue[13] = par[39];  customvalue[14] = par[18];
            // dir
            customvalue[15] = par[28]; customvalue[16] = par[29];
            // load
            customvalue[17] = par[30]; customvalue[18] = par[31];

            for(int i = 0; i < 19; i++) {
                printf("customvalue[%d] = %s = %lf\n", i, names[i], customvalue[i]);
            }
        }
    }
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) { 
        return;
    }
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        // Sigma
            customvalue[0] = par[6];  customvalue[1] = par[7];   customvalue[2] = par[8];
            // Strain
            customvalue[3] = par[9];  customvalue[4] = par[10];  customvalue[5] = par[11];
            // Beta
            customvalue[6] = par[12]; customvalue[7] = par[13];  customvalue[8] = par[14];
            // E
            customvalue[9] = par[15]; customvalue[10] = par[38];  customvalue[11] = par[16];
            // Alpha
            customvalue[12] = par[17]; customvalue[13] = par[39];  customvalue[14] = par[18];
            // dir
            customvalue[15] = par[28]; customvalue[16] = par[29];
            // load
            customvalue[17] = par[30]; customvalue[18] = par[31];
    } 
    // Mode to perform calculations at the end of the time series
    else if (mode == 3) {
        return;
    }
    else {
        print_debug("Custom Function using mode = %d, please use 0, 1, 2 or 3", mode);
    }
}

void customcalc_pendulum_EMEH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
	return;
}

void customcalc_linear_oscillator_gravity(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
	return;
}

void customcalc_chuas_circuit(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
	return;
}

void customcalc_multidirectional_hybrid_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
	return;
}

void customcalc_multidirectional_hybrid_EH_zero_pend_length(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
	return;
}

void customcalc_pendulum_EMEH_dimensional(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
	return;
}

void customcalc_tetrastable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode) {
	return;
}

/* Model for customcalc functions: 

void customcalc_name(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, double *customvalue, int mode)
    // Mode to define names to be printed in the output file
    if (mode == 0) {    
        char *names[] = {   "list", "of", "names"   };
        // Assign names to custom values
        assign_names(names, ncustomvalues, customnames);
    } 
    // Mode to perform calculations in steady state regime of the time series
    else if (mode == 1) {
        // Add calculations as
        customvalue[0] = operations;
        customvalue[1] = operations1;
    } 
    // Mode to perform calculations over the entire time series (transient + steady state)
    else if (mode == 2) {
        // Add calculations as
        customvalue[2] = operations2;
        customvalue[3] = operations3;
    }
    // Mode to perform calculations at the end of the time series    
    else if (mode == 3) {
        // Add calculations as
        customvalue[4] = operations4;
        customvalue[5] = operations5;
    }
    else {
        error(mode);
    }
}
*/

