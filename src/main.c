#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libs/interface.h"
#include "libs/odesystems.h"
#include "libs/customcalc.h"
#include "modules/time_series.h"
#include "modules/poinc_map.h"
#include "modules/lyap_exp_wolf.h"
#include "modules/ftime_series.h"
#include "modules/bifurcation.h"
#include "modules/fbifurcation.h"
#include "modules/dyndiag.h"
#include "modules/fdyndiag.h"
#include "modules/epbasin.h"
#include "modules/forcedbasin.h"
#include "modules/OS_time_series.h"
#include "modules/OS_ftime_series.h"
#include "modules/OS_bifurcation.h"
#include "modules/OS_fbifurcation.h"
#include "modules/OS_dyndiag.h"
#include "modules/OS_fdyndiag.h"
#include "modules/OS_fforcedbasin.h"
#include "modules/convergence_test.h"

void execute_OS_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), 
                        void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, size_t, double *, int),
                        char* outputname, char *funcname);
void execute_GNL_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname, char *funcname);

#define GNL_FUNC_1 lorenz
#define GNL_OUTPUTNAME_1 "lorenz"
#define GNL_FUNC_2 lotka_volterra_predator_prey
#define GNL_OUTPUTNAME_2 "lotka-volterra"
#define GNL_FUNC_3 halvorsen
#define GNL_OUTPUTNAME_3 "halvorsen"

#define OS_FUNC_1 duffing
#define OS_CUSTOM_1 customcalc
#define OS_OUTPUTNAME_1 "duffing"

#define OS_FUNC_2 duffing_2DoF                 
#define OS_CUSTOM_2 customcalc
#define OS_OUTPUTNAME_2 "duffing_2DoF"

#define OS_FUNC_3 vanderpol
#define OS_CUSTOM_3 customcalc
#define OS_OUTPUTNAME_3 "vanderpol"

#define OS_FUNC_4 pendulum
#define OS_CUSTOM_4 customcalc
#define OS_OUTPUTNAME_4 "pendulum"

#define OS_FUNC_5 falksma
#define OS_CUSTOM_5 customcalc
#define OS_OUTPUTNAME_5 "falk_sma"

#define OS_FUNC_6 linear_oscillator
#define OS_CUSTOM_6 customcalc
#define OS_OUTPUTNAME_6 "linear_oscillator"

#define OS_FUNC_7 duffing_vanderpol
#define OS_CUSTOM_7 customcalc
#define OS_OUTPUTNAME_7 "duffing_vanderpol"

#define OS_FUNC_8 bistable_EH               
#define OS_CUSTOM_8 customcalc_bistable_EH
#define OS_OUTPUTNAME_8 "bistable_EH"

#define OS_FUNC_9 tristable_EH
#define OS_CUSTOM_9 customcalc               
#define OS_OUTPUTNAME_9 "tristable_EH"

#define OS_FUNC_10 pend_oscillator_EH
#define OS_CUSTOM_10 customcalc_pend_oscillator_EH                 
#define OS_OUTPUTNAME_10 "pend_oscillator_EH"

#define OS_FUNC_11 pend_oscillator_wout_pend_EH
#define OS_CUSTOM_11 customcalc                 
#define OS_OUTPUTNAME_11 "pend_oscillator_EH(without_pend)"

#define OS_FUNC_12 duffing_2DoF_EH
#define OS_CUSTOM_12 customcalc                 
#define OS_OUTPUTNAME_12 "duffing_2DoF_EH"

#define OS_FUNC_13 linear_oscillator_2DoF
#define OS_CUSTOM_13 customcalc                 
#define OS_OUTPUTNAME_13 "lin_oscillator_2DoF"

#define MAX_NAMELENGTH 120
#define NUM_OF_GNL_SYSTEMS 3
#define NUM_OF_OS_SYSTEMS 13
#define NUM_OF_TOOLBOXES 2
#define NUM_OF_OS_MODULES 11
#define NUM_OF_GNL_MODULES 3 

char *OSsystemNames[NUM_OF_OS_SYSTEMS] = {  "Duffing Oscillator",
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
                                            "2 DoF Linear Oscillator" };
 
char *GNLsystemNames[NUM_OF_GNL_SYSTEMS] = {"Lorenz System",
                                            "Lotka-Volterra Predator-Prey Model",
                                            "Halvorsen System" };

char *toolboxesNames[NUM_OF_TOOLBOXES] = {  "General Nonlinear Dynamics Toolbox",
                                            "Harmonic Nonlinear Oscillators Toolbox" };
                            
char *OSmoduleNames[NUM_OF_OS_MODULES] = {  "Convergence Test",
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

char *GNLmoduleNames[NUM_OF_GNL_MODULES] = {"Convergence Test",
                                            "Time Series",
                                            "Lyapunov Exponents (Method from Wolf et al., 1985)"};

void call_GNL_system(unsigned int system, unsigned int module) {
        switch(system) {
        case 1:
            execute_GNL_modules(module, GNL_FUNC_1, GNL_OUTPUTNAME_1, GNLsystemNames[0]);
            break;
        case 2:
            execute_GNL_modules(module, GNL_FUNC_2, GNL_OUTPUTNAME_2, GNLsystemNames[1]);
            break;
        case 3:
            execute_GNL_modules(module, GNL_FUNC_3, GNL_OUTPUTNAME_3, GNLsystemNames[2]);
            break;
        default:
            printf("Invalid...\n");
            exit(0);
    }
}

void call_OS_system(unsigned int system, unsigned int module) {
    switch(system) {
        case 1:
            execute_OS_modules(module, OS_FUNC_1, OS_CUSTOM_1, OS_OUTPUTNAME_1, OSsystemNames[0]);
            break;
        case 2:
            execute_OS_modules(module, OS_FUNC_2, OS_CUSTOM_2, OS_OUTPUTNAME_2, OSsystemNames[1]);
            break;
        case 3:
            execute_OS_modules(module, OS_FUNC_3, OS_CUSTOM_3, OS_OUTPUTNAME_3, OSsystemNames[2]);
            break;
        case 4:
            execute_OS_modules(module, OS_FUNC_4, OS_CUSTOM_4, OS_OUTPUTNAME_4, OSsystemNames[3]);
            break;
        case 5:
            execute_OS_modules(module, OS_FUNC_5, OS_CUSTOM_5, OS_OUTPUTNAME_5, OSsystemNames[4]);
            break;
        case 6:
            execute_OS_modules(module, OS_FUNC_6, OS_CUSTOM_6, OS_OUTPUTNAME_6, OSsystemNames[5]);
            break;
        case 7:
            execute_OS_modules(module, OS_FUNC_7, OS_CUSTOM_7, OS_OUTPUTNAME_7, OSsystemNames[6]);
            break;
        case 8:
            execute_OS_modules(module, OS_FUNC_8, OS_CUSTOM_8, OS_OUTPUTNAME_8, OSsystemNames[7]);
            break;
        case 9:
            execute_OS_modules(module, OS_FUNC_9, OS_CUSTOM_9, OS_OUTPUTNAME_9, OSsystemNames[8]);
            break;
        case 10:
            execute_OS_modules(module, OS_FUNC_10, OS_CUSTOM_10, OS_OUTPUTNAME_10, OSsystemNames[9]);
            break;
        case 11:
            execute_OS_modules(module, OS_FUNC_11, OS_CUSTOM_11, OS_OUTPUTNAME_11, OSsystemNames[10]);
            break;
        case 12:
            execute_OS_modules(module, OS_FUNC_12, OS_CUSTOM_12, OS_OUTPUTNAME_12, OSsystemNames[11]);
            break;
        case 13:
        execute_OS_modules(module, OS_FUNC_13, OS_CUSTOM_13, OS_OUTPUTNAME_13, OSsystemNames[12]);
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
        identify_simulation(toolbox, &system, &module, toolboxesNames, OSsystemNames, OSmoduleNames, NUM_OF_OS_SYSTEMS, MAX_NAMELENGTH, NUM_OF_OS_MODULES);
        call_OS_system(system, module);
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

void execute_OS_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), 
                        void (*customfunc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, size_t, double *, int),
                        char* outputname, char *funcname) {
    switch (module) {
        case 1:
            convergence_test(funcname, outputname, edosys);
            break;
        case 2:
            OS_timeseries(funcname, outputname, edosys, customfunc);
            break;
        case 3:
            poincaremap(funcname, outputname, edosys);
            break;
        case 4:
            lyapunov_exp_wolf(funcname, outputname, edosys);
            break;
        case 5:
            OS_ftime_series(funcname, outputname, edosys, customfunc);
            break;
        case 6:
            OS_bifurcation(funcname, outputname, edosys, customfunc);
            break;
        case 7:
            OS_fbifurcation(funcname, outputname, edosys, customfunc);
            break;
        case 8:
            OS_dyndiag(funcname, outputname, edosys, customfunc);
            break;
        case 9:
            OS_fdyndiag(funcname, outputname, edosys, customfunc);
            break;
        case 10:
            epbasin(funcname, outputname, edosys);
            break;
        case 11:
            OS_fforcedbasin(funcname, outputname, edosys, customfunc);
            break;    
        default:
            printf("Invalid Module\n");
            exit(0);
    }
}

void execute_GNL_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname, char *funcname) {
    switch (module) {
        case 1:
            convergence_test(funcname, outputname, edosys);
            break;
        case 2:
            timeseries(funcname, outputname, edosys);
            break;
        case 3:
            lyapunov_exp_wolf(funcname, outputname, edosys);
            break;
        default:
            printf("Invalid Module\n");
            exit(0);
    }
}