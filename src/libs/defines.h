#ifndef DEFINES_H
#define DEFINES_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

// macros for compilers
#define _CRT_SECURE_NO_WARNINGS

// macros to handle file names
#define MAX_FILE_NUM 1000           // maximum number of repeatable files that the program can create in the directory
#define MAX_FILENAME_LEN 500        // Define maximum length of filename strings

// Macros related to printing information
#define MAX_PRINT_LEN 71            // Max length of the info printed on the screen and on info file
#define PERC_PRINT_NAME 0.6         // Percentage of space occuped by the name of the quantity printed
#define MAX_CCALC_NAME_LEN 31       // Maximum length of custom calculations names
#define PROGRESS_BAR_LEN 49         // Length of the progress bar
#define LF "%.10lf"                 // Significant digits for double numbers
#define SC "%.10e"                  // Signigicant digits for double numbers in scientific format 

// Math constants
//#define PI (4.0 * atan(1.0))
#define PI 3.141592653589793
#define TWOPI 2.0*PI

// Macro to handle variable names
#define getName(var) #var

// Define struct to handle if system has any angles
typedef struct {
    unsigned int n_angles;
    unsigned int *index;
} ang_info;

// Define struct to handle timeseries parameters (still unused)
typedef struct {
    unsigned int np;
    unsigned int ndiv;
    double tf;
    double h;
    int trans;
} ts;

// Define struct to handle initial conditions (still unused)
typedef struct {
    double *x;
    double t0;
} IC;

// Define struct to handle lyapunov exponents (still unused)
typedef struct {
    unsigned int ndim;
    double *lambda;
    double *s_lambda;
    double s_t0;
    double *LE;
    int N;
    double steadystateperc;
} lyap;

// Define struct to handle rms values of state variables (still unused)
typedef struct {
    double *xrms;
    double *overallxrms;
    int *rmsindex;
} rms;

// Define struct to handle determination of dynamical attractors (still unused)
typedef struct {
    int attractor;
    double **poinc;
    int numtol;
} attr;

// Define struct to handle mininum and maximum values (still unused)
typedef struct {
    double *xmin;
    double *xmax;
    double *overallxmin;
    double *overallxmax;
} minmax;

// Define struct to handle 1D diagrams params (still unused)
typedef struct {
    int index;
    double cpar_step;
    double *cpar_range;
    int mode;
} diag1D;

// Define struct to handle 2D diagrams params (still unused)
typedef struct {
    int index_X;
    double cpar_step_X;
    double *cpar_range_X;
    int index_Y;
    double cpar_step_Y;
    double *cpar_range_Y;
    int mode;
} diag2D;

// Define struct to handle custom calculations (still unused)
typedef struct {
    double value;
    char *name;
    unsigned int index;
    bool print;
    int position;
} cvar;

// Define struct to handle files (still unused)
typedef struct {
    FILE *results;
    FILE *poinc_results;
    FILE *info;
} files;

// main program macros 
#define MAX_NAMELENGTH 120
#define NUM_OF_GNL_SYSTEMS 4
#define NUM_OF_HOS_SYSTEMS 23
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
#define HOS_ANGLES_1 0

#define HOS_FUNC_2 duffing_2DoF                 
#define HOS_CUSTOM_2 customcalc
#define HOS_OUTPUTNAME_2 "duffing_2DoF"
#define HOS_DIM_2 4
#define HOS_NPAR_2 9
#define HOS_ANGLES_2 0

#define HOS_FUNC_3 vanderpol
#define HOS_CUSTOM_3 customcalc
#define HOS_OUTPUTNAME_3 "vanderpol"
#define HOS_DIM_3 2
#define HOS_NPAR_3 3
#define HOS_ANGLES_3 0

#define HOS_FUNC_4 pendulum
#define HOS_CUSTOM_4 customcalc
#define HOS_OUTPUTNAME_4 "pendulum"
#define HOS_DIM_4 2
#define HOS_NPAR_4 5
#define HOS_ANGLES_4 1
#define HOS_ANGINDEX0_4 0

#define HOS_FUNC_5 falksma
#define HOS_CUSTOM_5 customcalc
#define HOS_OUTPUTNAME_5 "falk_sma"
#define HOS_DIM_5 2
#define HOS_NPAR_5 6
#define HOS_ANGLES_5 0

#define HOS_FUNC_6 linear_oscillator
#define HOS_CUSTOM_6 customcalc
#define HOS_OUTPUTNAME_6 "linear_oscillator"
#define HOS_DIM_6 2
#define HOS_NPAR_6 4
#define HOS_ANGLES_6 0

#define HOS_FUNC_7 duffing_vanderpol
#define HOS_CUSTOM_7 customcalc
#define HOS_OUTPUTNAME_7 "duffing_vanderpol"
#define HOS_DIM_7 2
#define HOS_NPAR_7 6
#define HOS_ANGLES_7 0

#define HOS_FUNC_8 bistable_EH               
#define HOS_CUSTOM_8 customcalc_bistable_EH
#define HOS_OUTPUTNAME_8 "bistable_EH"
#define HOS_DIM_8 3
#define HOS_NPAR_8 8
#define HOS_ANGLES_8 0

#define HOS_FUNC_9 tristable_EH
#define HOS_CUSTOM_9 customcalc_tristable_EH               
#define HOS_OUTPUTNAME_9 "tristable_EH"
#define HOS_DIM_9 3
#define HOS_NPAR_9 9
#define HOS_ANGLES_9 0

#define HOS_FUNC_10 pend_oscillator_EH
#define HOS_CUSTOM_10 customcalc_pend_oscillator_EH                 
#define HOS_OUTPUTNAME_10 "pend_oscillator_EH"
#define HOS_DIM_10 8
#define HOS_NPAR_10 17
#define HOS_ANGLES_10 1
#define HOS_ANGINDEX0_10 4

#define HOS_FUNC_11 pend_oscillator_wout_pend_EH
#define HOS_CUSTOM_11 customcalc                 
#define HOS_OUTPUTNAME_11 "pend_oscillator_EH(without_pend)"
#define HOS_DIM_11 5
#define HOS_NPAR_11 10
#define HOS_ANGLES_11 0

#define HOS_FUNC_12 duffing_2DoF_EH
#define HOS_CUSTOM_12 customcalc_duffing_2DoF_EH                 
#define HOS_OUTPUTNAME_12 "duffing_2DoF_EH"
#define HOS_DIM_12 6
#define HOS_NPAR_12 16
#define HOS_ANGLES_12 0

#define HOS_FUNC_13 linear_oscillator_2DoF
#define HOS_CUSTOM_13 customcalc                 
#define HOS_OUTPUTNAME_13 "lin_oscillator_2DoF"
#define HOS_DIM_13 4
#define HOS_NPAR_13 6
#define HOS_ANGLES_13 0

#define HOS_FUNC_14 linear_2DoF_EH
#define HOS_CUSTOM_14 customcalc_linear_2DoF_EH                 
#define HOS_OUTPUTNAME_14 "lin_2DoF_EH"
#define HOS_DIM_14 6
#define HOS_NPAR_14 12
#define HOS_ANGLES_14 0

#define HOS_FUNC_15 adeodato_sma_oscillator
#define HOS_CUSTOM_15 customcalc_adeodato_sma_oscillator                 
#define HOS_OUTPUTNAME_15 "adeodato_sma_oscillator"
#define HOS_DIM_15 2 
#define HOS_NPAR_15 40
#define HOS_ANGLES_15 0

#define HOS_FUNC_16 linear_EMEH
#define HOS_CUSTOM_16 customcalc                 
#define HOS_OUTPUTNAME_16 "lin_EM_EH"
#define HOS_DIM_16 3
#define HOS_NPAR_16 9
#define HOS_ANGLES_16 0

#define HOS_FUNC_17 pendulum_EMEH
#define HOS_CUSTOM_17 customcalc_pendulum_EMEH
#define HOS_OUTPUTNAME_17 "pendulum_EMEH"
#define HOS_DIM_17 3
#define HOS_NPAR_17 7
#define HOS_ANGLES_17 1
#define HOS_ANGINDEX0_17 0

#define HOS_FUNC_18 linear_oscillator_gravity
#define HOS_CUSTOM_18 customcalc_linear_oscillator_gravity
#define HOS_OUTPUTNAME_18 "lin_oscillator_gravity"
#define HOS_DIM_18 2
#define HOS_NPAR_18 6
#define HOS_ANGLES_18 0

#define HOS_FUNC_19 multidirectional_hybrid_EH
#define HOS_CUSTOM_19 customcalc_multidirectional_hybrid_EH                 
#define HOS_OUTPUTNAME_19 "multidirect_hybrid_EH"
#define HOS_DIM_19 8
#define HOS_NPAR_19 16
#define HOS_ANGLES_19 1
#define HOS_ANGINDEX0_19 4

#define HOS_FUNC_20 multidirectional_hybrid_EH_zero_pend_length
#define HOS_CUSTOM_20 customcalc_multidirectional_hybrid_EH_zero_pend_length                 
#define HOS_OUTPUTNAME_20 "multidirect_hybrid_EH_zero_len_pend"
#define HOS_DIM_20 5
#define HOS_NPAR_20 10
#define HOS_ANGLES_20 0

#define HOS_FUNC_21 pendulum_EMEH_dimensional
#define HOS_CUSTOM_21 customcalc_pendulum_EMEH_dimensional
#define HOS_OUTPUTNAME_21 "pendulum_EMEH_dimensional"
#define HOS_DIM_21 3
#define HOS_NPAR_21 7
#define HOS_ANGLES_21 1
#define HOS_ANGINDEX0_21 0

#define HOS_FUNC_22 tetrastable_EH
#define HOS_CUSTOM_22 customcalc_tetrastable_EH               
#define HOS_OUTPUTNAME_22 "tetrastable_EH"
#define HOS_DIM_22 3
#define HOS_NPAR_22 10
#define HOS_ANGLES_22 0

#define HOS_FUNC_23 multidirectional_hybrid_EH_coupling_ratio
#define HOS_CUSTOM_23 customcalc_multidirectional_hybrid_EH_coupling_ratio                 
#define HOS_OUTPUTNAME_23 "multidirect_hybrid_EH_CR"
#define HOS_DIM_23 8
#define HOS_NPAR_23 15
#define HOS_ANGLES_23 1
#define HOS_ANGINDEX0_23 4

#endif
