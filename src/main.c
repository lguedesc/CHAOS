#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libs/interface.h"
#include "libs/odesystems.h"
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
#include "modules/EH_time_series.h"
#include "modules/EH_ftime_series.h"
#include "modules/EH_bifurcation.h"
#include "modules/EH_fbifurcation.h"
#include "modules/EH_dyndiag.h"
#include "modules/EH_fdyndiag.h"
#include "modules/EH_forcedbasin.h"
#include "modules/convergence_test.h"

void execute_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname, char *funcname);
void execute_EH_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname, char *funcname);

#define FUNC_1 duffing
#define OUTPUTNAME_1 "duffing"
#define FUNC_2 duffing_2DoF                 
#define OUTPUTNAME_2 "duffing_2DoF"
#define FUNC_3 vanderpol
#define OUTPUTNAME_3 "vanderpol"
#define FUNC_4 pendulum
#define OUTPUTNAME_4 "pendulum"
#define FUNC_5 falksma
#define OUTPUTNAME_5 "falk_sma"
#define FUNC_6 linear_oscillator
#define OUTPUTNAME_6 "linear_oscillator"
#define FUNC_7 duffing_vanderpol
#define OUTPUTNAME_7 "duffing_vanderpol"

#define EH_FUNC_1 bistable_EH               
#define EH_OUTPUTNAME_1 "bistable_EH"
#define EH_FUNC_2 tristable_EH               
#define EH_OUTPUTNAME_2 "tristable_EH"
#define EH_FUNC_3 pend_oscillator_EH               
#define EH_OUTPUTNAME_3 "pend_oscillator_EH"
#define EH_FUNC_4 pend_oscillator_wout_pend_EH               
#define EH_OUTPUTNAME_4 "pend_oscillator_EH(without_pend)"

#define MAX_NAMELENGTH 120
#define NUM_OF_SYSTEMS 7
#define NUM_OF_EH_SYSTEMS 4
#define NUM_OF_TOOLBOXES 2
#define NUM_OF_MODULES 11 

char *systemNames[NUM_OF_SYSTEMS] = {"Duffing Oscillator",
                                     "2 DoF Duffing Oscillator",
                                     "Van Der Pol Oscillator",
                                     "Simple Pendulum",
                                     "Falk Shape Memory Alloy Oscillator",
                                     "Linear Oscillator",
                                     "Duffing-Van Der Pol Oscillator"};
 
char *EHsystemNames[NUM_OF_EH_SYSTEMS] = {"Polynomial Bistable Energy Harvester",
                                          "Polynomial Tristable Energy Harvester",
                                          "Pendulum-Oscillator Energy Harvester",
                                          "Pendulum-Oscillator Energy Harvester (Without Pendulum)"};

char *toolboxesNames[NUM_OF_TOOLBOXES] = {"Nonlinear Dynamics Toolbox",
                                          "Energy Harvesting Toolbox"};
                            
char *moduleNames[NUM_OF_MODULES] = {"Convergence Test",
                                     "Time Series",
                                     "Poincare Map",
                                     "Lyapunov Exponents (Method from Wolf et al., 1985)", 
                                     "Full Time Series (Integrator + Poincare Map + Lyapunov Exponents)",
                                     "Bifurcation Diagram",
                                     "Full Bifurcation Diagram (With Lyapunov)",
                                     "Dynamical Diagram",
                                     "Full Dynamical Diagram (With Lyapunov)",
                                     "Basin of Attraction (Fixed Points)",
                                     "Basin of Attraction (Forced)"};

void call_system(unsigned int system, unsigned int module) {
    switch(system) {
        case 1:
            execute_modules(module, FUNC_1, OUTPUTNAME_1, systemNames[0]);
            break;
        case 2:
            execute_modules(module, FUNC_2, OUTPUTNAME_2, systemNames[1]);
            break;
        case 3:
            execute_modules(module, FUNC_3, OUTPUTNAME_3, systemNames[2]);
            break;
        case 4:
            execute_modules(module, FUNC_4, OUTPUTNAME_4, systemNames[3]);
            break;
        case 5:
            execute_modules(module, FUNC_5, OUTPUTNAME_5, systemNames[4]);
            break;
        case 6:
            execute_modules(module, FUNC_6, OUTPUTNAME_6, systemNames[5]);
            break;
        case 7:
            execute_modules(module, FUNC_7, OUTPUTNAME_7, systemNames[6]);
            break;  
        default:
            printf("Invalid...\n");
            exit(0);
    }
}

void call_EH_system(unsigned int system, unsigned int module) {
    switch(system) {
        case 1:
            execute_EH_modules(module, EH_FUNC_1, EH_OUTPUTNAME_1, EHsystemNames[0]);
            break;
        case 2:
            execute_EH_modules(module, EH_FUNC_2, EH_OUTPUTNAME_2, EHsystemNames[1]);
            break;
        case 3:
            execute_EH_modules(module, EH_FUNC_3, EH_OUTPUTNAME_3, EHsystemNames[2]);
            break;
        case 4:
            execute_EH_modules(module, EH_FUNC_4, EH_OUTPUTNAME_4, EHsystemNames[3]);
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
        identify_simulation(toolbox, &system, &module, toolboxesNames, systemNames, moduleNames, NUM_OF_SYSTEMS, MAX_NAMELENGTH, NUM_OF_MODULES);
        call_system(system, module);
    }
    else if (toolbox == 2) {
        identify_simulation(toolbox, &system, &module, toolboxesNames, EHsystemNames, moduleNames, NUM_OF_EH_SYSTEMS, MAX_NAMELENGTH, NUM_OF_MODULES);
        call_EH_system(system, module);
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

void execute_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname, char *funcname) {
    switch (module) {
        case 1:
            convergence_test(funcname, outputname, edosys);
            break;
        case 2:
            timeseries(funcname, outputname, edosys);
            break;
        case 3:
            poincaremap(funcname, outputname, edosys);
            break;
        case 4:
            lyapunov_exp_wolf(funcname, outputname, edosys);
            break;
        case 5:
            ftime_series(funcname, outputname, edosys);
            break;
        case 6:
            bifurcation(funcname, outputname, edosys);
            break;
        case 7:
            fbifurcation(funcname, outputname, edosys);
            break;
        case 8:
            dyndiag(funcname, outputname, edosys);
            break;
        case 9:
            fdyndiag(funcname, outputname, edosys);
            break;
        case 10:
            epbasin(funcname, outputname, edosys);
            break;
        case 11:
            forcedbasin(funcname, outputname, edosys);
            break;    
        default:
            printf("Invalid Module\n");
            exit(0);
    }
}

void execute_EH_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname, char *funcname) {
    switch (module) {
        case 1:
            convergence_test(funcname, outputname, edosys);
            break;
        case 2:
            EH_timeseries(funcname, outputname, edosys);
            break;
        case 3:
            poincaremap(funcname, outputname, edosys);
            break;
        case 4:
            lyapunov_exp_wolf(funcname, outputname, edosys);
            break;
        case 5:
            EH_ftime_series(funcname, outputname, edosys);
            break;
        case 6:
            EH_bifurcation(funcname, outputname, edosys);
            break;
        case 7:
            EH_fbifurcation(funcname, outputname, edosys);
            break;
        case 8:
            EH_dyndiag(funcname, outputname, edosys);
            break;
        case 9:
            EH_fdyndiag(funcname, outputname, edosys);
            break;
        case 10:
            epbasin(funcname, outputname, edosys);
            break;
        case 11:
            EH_forcedbasin(funcname, outputname, edosys);
            break;    
        default:
            printf("Invalid Module\n");
            exit(0);
    }
}