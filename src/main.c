#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libs/interface.h"
#include "libs/edosystems.h"
#include "modules/time_series.h"
#include "modules/poinc_map.h"
#include "modules/lyap_exp_wolf.h"
#include "modules/ftime_series.h"
#include <crtdbg.h>

char **create_memory_to_store_names(size_t rows, size_t cols);
void execute_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname);

#define FUNC_1 duffing
#define OUTPUTNAME_1 "duffing"
#define FUNC_2 duffing_2DoF                 
#define OUTPUTNAME_2 "duffing_2DoF"

#define EH_FUNC_1 bistable_EH               
#define EH_OUTPUTNAME_1 "bistable_EH"

#define MAX_NAMELENGTH 100
#define NUM_OF_SYSTEMS 2
#define NUM_OF_EH_SYSTEMS 1
#define NUM_OF_TOOLBOXES 2
#define NUM_OF_MODULES 9

void call_system(unsigned int system, unsigned int module) {
    switch(system) {
        case 1:
            execute_modules(module, FUNC_1, OUTPUTNAME_1);
            break;
        case 2:
            execute_modules(module, FUNC_2, OUTPUTNAME_2);
            break;
        default:
            printf("Invalid...\n");
            exit(0);
    }
}

void call_EH_system(unsigned int system, unsigned int module) {
    switch(system) {
        case 1:
            execute_modules(module, EH_FUNC_1, EH_OUTPUTNAME_1);
            break;
        default:
            printf("Invalid...\n");
            exit(0);
    }
}

void assign_system_names(char **systemnames, char **EHsystemnames) {
    // Nonlinear Dynamics Toolbox System Names
    strcpy(systemnames[0], "Duffing Oscillator");
    strcpy(systemnames[1], "2 DoF Duffing Oscillator");
    // Energy Harvesting Toolbox System Names
    strcpy(EHsystemnames[0], "Bistable Energy Harvester");
}

void assign_toolbox_names(char **toolboxnames) {
    // Toolbox Names
    strcpy(toolboxnames[0], "Nonlinear Dynamics Toolbox");
    strcpy(toolboxnames[1], "Energy Harvesting Toolbox");
}

void assign_module_names(char **modulenames) {
    // Module Names
    strcpy(modulenames[0], "Time Series");
    strcpy(modulenames[1], "Poincare Map");
    strcpy(modulenames[2], "Lyapunov Exponents (Method from Wolf et al., 1985)");
    strcpy(modulenames[3], "Full Time Series (Integrator + Poincare Map + Lyapunov Exponents)");
    strcpy(modulenames[4], "Bifurcation Diagram");
    strcpy(modulenames[5], "Full Bifurcation Diagram (Automatic Identification of Attractors)");
    strcpy(modulenames[6], "Dynamical Diagram");
    strcpy(modulenames[7], "Basin of Attraction (Fixed Points)");
    strcpy(modulenames[8], "Basin of Attraction (Forced)");
}


int main (void) {

    char **systemNames = create_memory_to_store_names(NUM_OF_SYSTEMS, MAX_NAMELENGTH);
    char **EHsystemNames = create_memory_to_store_names(NUM_OF_EH_SYSTEMS, MAX_NAMELENGTH);
    char **toolboxesNames = create_memory_to_store_names(NUM_OF_TOOLBOXES, MAX_NAMELENGTH);
    char **moduleNames = create_memory_to_store_names(NUM_OF_MODULES, MAX_NAMELENGTH);
    assign_system_names(systemNames, EHsystemNames);
    assign_toolbox_names(toolboxesNames);
    assign_module_names(moduleNames);
    
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
    
    // Free Memory
    free(systemNames);
    free(EHsystemNames);
    free(toolboxesNames);
    free(moduleNames);
    
    // End of execution
    end_of_execution(MAX_NAMELENGTH);
}

char **create_memory_to_store_names(size_t rows, size_t cols) {
    char **arr = malloc(rows * sizeof **arr);
    for (int i = 0; i < rows; i++) {
        arr[i] = malloc(cols * sizeof **arr);
    }
    return arr;
}

void execute_modules(unsigned int module, void (*edosys)(int, double *, double, double *, double *), char* outputname) {
    switch (module) {
        case 1:
            timeseries(outputname, edosys);
            break;
        case 2:
            poincaremap(outputname, edosys);
            break;
        case 3:
            lyapunov_exp_wolf(outputname, edosys);
            break;
        case 4:
            ftime_series(outputname, edosys);
            break;
        case 5:
            // Call solution of the type solution(edosys, outputname);
            break;
        case 6:
            // Call solution of the type solution(edosys, outputname);
            break;
        case 7:
            // Call solution of the type solution(edosys, outputname);
            break;
        case 8:
            // Call solution of the type solution(edosys, outputname);
            break;
        case 9:
            // Call solution of the type solution(edosys, outputname);
            break;    
        default:
            printf("Invalid Module\n");
            exit(0);
    }
}

