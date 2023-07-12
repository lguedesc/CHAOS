#ifndef DEFINES_H
#define DEFINES_H

// macros for compilers
#define _CRT_SECURE_NO_WARNINGS

// macros to handle file names
#define MAX_FILE_NUM 1000           // maximum number of repeatable files that the program can create in the directory
#define MAX_FILENAME_LEN 500        // Define maximum length of filename strings

// Macros related to printing information
#define MAX_PRINT_LEN 71            // Max length of the info printed on the screen and on info file
#define PERC_PRINT_NAME 0.6         // Percentage of space occuped by the name of the quantity printed

// Math constants
#define PI (4 * atan(1))

// main program macros 
#define MAX_NAMELENGTH 120
#define NUM_OF_GNL_SYSTEMS 4
#define NUM_OF_HOS_SYSTEMS 16
#define NUM_OF_TOOLBOXES 2
#define NUM_OF_HOS_MODULES 11
#define NUM_OF_GNL_MODULES 3

// macros for general nonlinear (GNL) ode system identification
#define GNL_FUNC_1 lorenz
#define GNL_CUSTOM_1 customcalc_lorenz
#define GNL_OUTPUTNAME_1 "lorenz"
#define GNL_DIM_1 3
#define GNL_NPAR_1 4  // 3 parameters + omega for determine stepsize

#define GNL_FUNC_2 lotka_volterra_predator_prey
#define GNL_CUSTOM_2 customcalc_lotka_volterra_predator_prey
#define GNL_OUTPUTNAME_2 "lotka-volterra"
#define GNL_DIM_2 2
#define GNL_NPAR_2 5  // 4 parameters + omega for determine stepsize

#define GNL_FUNC_3 halvorsen
#define GNL_CUSTOM_3 customcalc_halvorsen
#define GNL_OUTPUTNAME_3 "halvorsen"
#define GNL_DIM_3 3
#define GNL_NPAR_3 2 // 1 parameter + omega for determine stepsize

#define GNL_FUNC_4 chuas_circuit
#define GNL_CUSTOM_4 customcalc_chuas_circuit
#define GNL_OUTPUTNAME_4 "chuas_circuit"
#define GNL_DIM_4 3
#define GNL_NPAR_4 5 // 4 parameters + omega for determine stepsize

// macros for general nonlinear (GNL) ode system identification
#define HOS_FUNC_1 duffing
#define HOS_CUSTOM_1 customcalc
#define HOS_OUTPUTNAME_1 "duffing"
#define HOS_DIM_1 2
#define HOS_NPAR_1 5
#define HOS_ANGLES_1 false

#define HOS_FUNC_2 duffing_2DoF                 
#define HOS_CUSTOM_2 customcalc
#define HOS_OUTPUTNAME_2 "duffing_2DoF"
#define HOS_DIM_2 4
#define HOS_NPAR_2 9
#define HOS_ANGLES_2 false

#define HOS_FUNC_3 vanderpol
#define HOS_CUSTOM_3 customcalc
#define HOS_OUTPUTNAME_3 "vanderpol"
#define HOS_DIM_3 2
#define HOS_NPAR_3 3
#define HOS_ANGLES_3 false

#define HOS_FUNC_4 pendulum
#define HOS_CUSTOM_4 customcalc
#define HOS_OUTPUTNAME_4 "pendulum"
#define HOS_DIM_4 2
#define HOS_NPAR_4 5
#define HOS_ANGLES_4 true

#define HOS_FUNC_5 falksma
#define HOS_CUSTOM_5 customcalc
#define HOS_OUTPUTNAME_5 "falk_sma"
#define HOS_DIM_5 2
#define HOS_NPAR_5 6
#define HOS_ANGLES_5 false

#define HOS_FUNC_6 linear_oscillator
#define HOS_CUSTOM_6 customcalc
#define HOS_OUTPUTNAME_6 "linear_oscillator"
#define HOS_DIM_6 2
#define HOS_NPAR_6 4
#define HOS_ANGLES_6 false

#define HOS_FUNC_7 duffing_vanderpol
#define HOS_CUSTOM_7 customcalc
#define HOS_OUTPUTNAME_7 "duffing_vanderpol"
#define HOS_DIM_7 2
#define HOS_NPAR_7 6
#define HOS_ANGLES_7 false

#define HOS_FUNC_8 bistable_EH               
#define HOS_CUSTOM_8 customcalc_bistable_EH
#define HOS_OUTPUTNAME_8 "bistable_EH"
#define HOS_DIM_8 3
#define HOS_NPAR_8 8
#define HOS_ANGLES_8 false

#define HOS_FUNC_9 tristable_EH
#define HOS_CUSTOM_9 customcalc               
#define HOS_OUTPUTNAME_9 "tristable_EH"
#define HOS_DIM_9 3
#define HOS_NPAR_9 9
#define HOS_ANGLES_9 false

#define HOS_FUNC_10 pend_oscillator_EH
#define HOS_CUSTOM_10 customcalc_pend_oscillator_EH                 
#define HOS_OUTPUTNAME_10 "pend_oscillator_EH"
#define HOS_DIM_10 8
#define HOS_NPAR_10 17
#define HOS_ANGLES_10 true

#define HOS_FUNC_11 pend_oscillator_wout_pend_EH
#define HOS_CUSTOM_11 customcalc                 
#define HOS_OUTPUTNAME_11 "pend_oscillator_EH(without_pend)"
#define HOS_DIM_11 5
#define HOS_NPAR_11 10
#define HOS_ANGLES_11 false

#define HOS_FUNC_12 duffing_2DoF_EH
#define HOS_CUSTOM_12 customcalc_duffing_2DoF_EH                 
#define HOS_OUTPUTNAME_12 "duffing_2DoF_EH"
#define HOS_DIM_12 6
#define HOS_NPAR_12 16
#define HOS_ANGLES_12 false

#define HOS_FUNC_13 linear_oscillator_2DoF
#define HOS_CUSTOM_13 customcalc                 
#define HOS_OUTPUTNAME_13 "lin_oscillator_2DoF"
#define HOS_DIM_13 4
#define HOS_NPAR_13 6
#define HOS_ANGLES_13 false

#define HOS_FUNC_14 linear_2DoF_EH
#define HOS_CUSTOM_14 customcalc_linear_2DoF_EH                 
#define HOS_OUTPUTNAME_14 "lin_2DoF_EH"
#define HOS_DIM_14 6
#define HOS_NPAR_14 12
#define HOS_ANGLES_14 false

#define HOS_FUNC_15 adeodato_sma_oscillator
#define HOS_CUSTOM_15 customcalc_adeodato_sma_oscillator                 
#define HOS_OUTPUTNAME_15 "adeodato_sma_oscillator"
#define HOS_DIM_15 2 
#define HOS_NPAR_15 40
#define HOS_ANGLES_15 false

#define HOS_FUNC_16 linear_EMEH
#define HOS_CUSTOM_16 customcalc                 
#define HOS_OUTPUTNAME_16 "lin_EM_EH"
#define HOS_DIM_16 3
#define HOS_NPAR_16 9
#define HOS_ANGLES_16 false

#endif