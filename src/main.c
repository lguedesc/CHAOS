#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "libs/interface.h"
#include "libs/odesystems.h"
#include "libs/customcalc.h"
#include "modules/time_series.h"
#include "modules/poinc_map.h"
#include "modules/lyap_exp_wolf.h"
#include "modules/epbasin.h"
#include "modules/HOS_timeseries.h"
#include "modules/HOS_ftime_series.h"
#include "modules/HOS_bifurcation.h"
#include "modules/HOS_fbifurcation.h"
#include "modules/HOS_dyndiag.h"
#include "modules/HOS_fdyndiag.h"
#include "modules/HOS_fforcedbasin.h"
#include "modules/convergence_test.h"
#include "libs/defines.h"

void execute_HOS_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), 
                        void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, size_t, double *, int),
                        char* outputname, char *funcname, unsigned int dim, unsigned int npar, bool angles);
void execute_GNL_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname,
                         char *funcname, unsigned int dim, unsigned int npar);

char *HOSsystemNames[NUM_OF_HOS_SYSTEMS] = {"Duffing Oscillator",
                                            "2 DoF Duffing Oscillator",
                                            "Van Der Pol Oscillator",
                                            "Simple Pendulum",
                                            "Falk Shape Memory Alloy Oscillator",
                                            "Linear Oscillator",
                                            "Duffing-Van Der Pol Oscillator",
                                            "Polynomial Bistable Energy Harvester",
                                            "Polynomial Tristable Energy Harvester",
                                            "Pendulum-Oscillator Energy Harvester",
                                            "Pendulum-Oscillator Energy Harvester (Without Pendulum)",
                                            "2 DoF Duffing-Type Energy Harvester",
                                            "2 DoF Linear Oscillator",
                                            "2 DoF Linear Energy Harvester",
                                            "Adeodato SMA Oscillator (In Development)",
                                            "Linear Electromagnetic Energy Harvester" };
 
char *GNLsystemNames[NUM_OF_GNL_SYSTEMS] = {"Lorenz System",
                                            "Lotka-Volterra Predator-Prey Model",
                                            "Halvorsen System",
											"Chua's Circuit" };

char *toolboxesNames[NUM_OF_TOOLBOXES] = {  "General Nonlinear Dynamics Toolbox",
                                            "Harmonic Nonlinear Oscillators Toolbox" };
                            
char *HOSmoduleNames[NUM_OF_HOS_MODULES] = {  "Convergence Test (In Development)",
                                            "Time Series",
                                            "Poincare Map",
                                            "Lyapunov Exponents (Method from Wolf et al., 1985)", 
                                            "Full Time Series (Integrator + Poincare Map + Lyapunov Exponents)",
                                            "Bifurcation Diagram",
                                            "Full Bifurcation Diagram (With Lyapunov)",
                                            "Dynamical Diagram",
                                            "Full Dynamical Diagram (With Lyapunov)",
                                            "Basin of Attraction (Fixed Points)",
                                            "Full Basin of Attraction (Forced)" };

char *GNLmoduleNames[NUM_OF_GNL_MODULES] = { "Convergence Test (In Development)",
                                             "Time Series",
                                             "Lyapunov Exponents (Method from Wolf et al., 1985)"};

void call_GNL_system(unsigned int system, unsigned int module) {
        switch(system) {
        case 1:
            execute_GNL_modules(module, GNL_FUNC_1, GNL_OUTPUTNAME_1, GNLsystemNames[0], GNL_DIM_1, GNL_NPAR_1);
            break;
        case 2:
            execute_GNL_modules(module, GNL_FUNC_2, GNL_OUTPUTNAME_2, GNLsystemNames[1], GNL_DIM_2, GNL_NPAR_2);
            break;
        case 3:
            execute_GNL_modules(module, GNL_FUNC_3, GNL_OUTPUTNAME_3, GNLsystemNames[2], GNL_DIM_3, GNL_NPAR_3);
            break;
		case 4:
			execute_GNL_modules(module, GNL_FUNC_4, GNL_OUTPUTNAME_4, GNLsystemNames[3], GNL_DIM_4, GNL_NPAR_4);
			break;

        default:
            printf("Invalid...\n");
            exit(0);
    }
}

void call_HOS_system(unsigned int system, unsigned int module) {
    switch(system) {
        case 1:
            execute_HOS_modules(module, HOS_FUNC_1, HOS_CUSTOM_1, HOS_OUTPUTNAME_1, HOSsystemNames[0], HOS_DIM_1, HOS_NPAR_1, HOS_ANGLES_1);
            break;
        case 2:
            execute_HOS_modules(module, HOS_FUNC_2, HOS_CUSTOM_2, HOS_OUTPUTNAME_2, HOSsystemNames[1], HOS_DIM_2, HOS_NPAR_2, HOS_ANGLES_2);
            break;
        case 3:
            execute_HOS_modules(module, HOS_FUNC_3, HOS_CUSTOM_3, HOS_OUTPUTNAME_3, HOSsystemNames[2], HOS_DIM_3, HOS_NPAR_3, HOS_ANGLES_3);
            break;
        case 4:
            execute_HOS_modules(module, HOS_FUNC_4, HOS_CUSTOM_4, HOS_OUTPUTNAME_4, HOSsystemNames[3], HOS_DIM_4, HOS_NPAR_4, HOS_ANGLES_4);
            break;
        case 5:
            execute_HOS_modules(module, HOS_FUNC_5, HOS_CUSTOM_5, HOS_OUTPUTNAME_5, HOSsystemNames[4], HOS_DIM_5, HOS_NPAR_5, HOS_ANGLES_5);
            break;
        case 6:
            execute_HOS_modules(module, HOS_FUNC_6, HOS_CUSTOM_6, HOS_OUTPUTNAME_6, HOSsystemNames[5], HOS_DIM_6, HOS_NPAR_6, HOS_ANGLES_6);
            break;
        case 7:
            execute_HOS_modules(module, HOS_FUNC_7, HOS_CUSTOM_7, HOS_OUTPUTNAME_7, HOSsystemNames[6], HOS_DIM_7, HOS_NPAR_7, HOS_ANGLES_7);
            break;
        case 8:
            execute_HOS_modules(module, HOS_FUNC_8, HOS_CUSTOM_8, HOS_OUTPUTNAME_8, HOSsystemNames[7], HOS_DIM_8, HOS_NPAR_8, HOS_ANGLES_8);
            break;
        case 9:
            execute_HOS_modules(module, HOS_FUNC_9, HOS_CUSTOM_9, HOS_OUTPUTNAME_9, HOSsystemNames[8], HOS_DIM_9, HOS_NPAR_9, HOS_ANGLES_9);
            break;
        case 10:
            execute_HOS_modules(module, HOS_FUNC_10, HOS_CUSTOM_10, HOS_OUTPUTNAME_10, HOSsystemNames[9], HOS_DIM_10, HOS_NPAR_10, HOS_ANGLES_10);
            break;
        case 11:
            execute_HOS_modules(module, HOS_FUNC_11, HOS_CUSTOM_11, HOS_OUTPUTNAME_11, HOSsystemNames[10], HOS_DIM_11, HOS_NPAR_11, HOS_ANGLES_11);
            break;
        case 12:
            execute_HOS_modules(module, HOS_FUNC_12, HOS_CUSTOM_12, HOS_OUTPUTNAME_12, HOSsystemNames[11], HOS_DIM_12, HOS_NPAR_12, HOS_ANGLES_12);
            break;
        case 13:
            execute_HOS_modules(module, HOS_FUNC_13, HOS_CUSTOM_13, HOS_OUTPUTNAME_13, HOSsystemNames[12], HOS_DIM_13, HOS_NPAR_13, HOS_ANGLES_13);
            break;
        case 14:
            execute_HOS_modules(module, HOS_FUNC_14, HOS_CUSTOM_14, HOS_OUTPUTNAME_14, HOSsystemNames[13], HOS_DIM_14, HOS_NPAR_14, HOS_ANGLES_14);
            break;
        case 15:
            execute_HOS_modules(module, HOS_FUNC_15, HOS_CUSTOM_15, HOS_OUTPUTNAME_15, HOSsystemNames[14], HOS_DIM_15, HOS_NPAR_15, HOS_ANGLES_15);
            break;
        case 16:
            execute_HOS_modules(module, HOS_FUNC_16, HOS_CUSTOM_16, HOS_OUTPUTNAME_16, HOSsystemNames[15], HOS_DIM_16, HOS_NPAR_16, HOS_ANGLES_16);
            break;
        default:
            printf("Invalid...\n");
            exit(0);
    }
}

int main (void) {

    // Define Control Variables
    unsigned int toolbox;
    unsigned int system;
    unsigned int module;
    
    welcome_header(MAX_NAMELENGTH);
    toolbox = choose_option(toolboxesNames, NUM_OF_TOOLBOXES, MAX_NAMELENGTH, "Toolbox");
    if (toolbox == 1) {
        identify_simulation(toolbox, &system, &module, toolboxesNames, GNLsystemNames, GNLmoduleNames, NUM_OF_GNL_SYSTEMS, MAX_NAMELENGTH, NUM_OF_GNL_MODULES);
        call_GNL_system(system, module);
    }
    else if (toolbox == 2) {
        identify_simulation(toolbox, &system, &module, toolboxesNames, HOSsystemNames, HOSmoduleNames, NUM_OF_HOS_SYSTEMS, MAX_NAMELENGTH, NUM_OF_HOS_MODULES);
        call_HOS_system(system, module);
    }
    else if (toolbox == 0) {
        clear_screen();
        exit(0);
    }
    else {
        clear_screen();
        welcome_header(MAX_NAMELENGTH);
        invalid_option(toolbox, "Toolbox", MAX_NAMELENGTH);
        partition(1, MAX_NAMELENGTH);
    }
    
    // End of execution
    end_of_execution(MAX_NAMELENGTH);
}

void execute_HOS_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), 
                        void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, size_t, double *, int),
                        char* outputname, char *funcname, unsigned int dim, unsigned int npar, bool angles) {
    switch (module) {
        case 1:
            convergence_test(funcname, dim, npar, outputname, edosys);
            break;
        case 2:
            HOS_timeseries(funcname, dim, npar, outputname, edosys, customfunc);
            break;
        case 3:
            poincaremap(funcname, dim, npar, outputname, edosys);  
            break;
        case 4:
            lyapunov_exp_wolf(funcname, dim, npar, outputname, edosys);
            break;
        case 5:
            HOS_ftime_series(funcname, dim, npar, outputname, edosys, customfunc);
            break;
        case 6:
            HOS_bifurcation(funcname, dim, npar, outputname, edosys, customfunc);
            break;
        case 7:
            HOS_fbifurcation(funcname, dim, npar, outputname, edosys, customfunc);
            break;
        case 8:
            HOS_dyndiag(funcname, dim, npar, outputname, edosys, customfunc);
            break;
        case 9:
            HOS_fdyndiag(funcname, dim, npar, outputname, edosys, customfunc);
            break;
        case 10:
            epbasin(funcname, dim, npar, outputname, edosys);
            break;
        case 11:
            HOS_fforcedbasin(funcname, dim, npar, outputname, edosys, customfunc);
            break;    
        default:
            printf("Invalid Module\n");
            exit(0);
    }
}

void execute_GNL_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname,
                         char *funcname, unsigned int dim, unsigned int npar) {
    switch (module) {
        case 1:
            convergence_test(funcname, dim, npar, outputname, edosys);
            break;
        case 2:
            timeseries(funcname, dim, npar, outputname, edosys);
            break;
        case 3:
            lyapunov_exp_wolf(funcname, dim, npar, outputname, edosys);
            break;
        default:
            printf("Invalid Module\n");
            exit(0);
    }
}