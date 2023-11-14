#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "msg.h"
#include "defines.h"

// Main Interface 
void clear_screen(void) {
    #ifdef _WIN32
        system("cls");
    #else
        system("clear");
    #endif
}

void partition(int mode, size_t maxlength) {
    if (mode == 1) {
        printf(" ");
        for (size_t i = 0; i < maxlength-2; i++) {
            printf("=");
        }
        printf(" ");
    }
    else if (mode == 2) {
        printf(" ");
        for (size_t i = 0; i < maxlength-2; i++) {
            printf("-");
        }
        printf(" ");
    }
    else {
        printf("DEBUG WARNING: Please select 1 for main partition or 2 for secondary partition for the interface programming");
        exit(0);
    }
    printf("\n");
}

void welcome_header(size_t maxlength) {
    //clear_screen();
    char *message = "Welcome to CHAOS, a Numerical Package for Nonlinear Dynamical Systems";
    partition(1, maxlength);
    int padlen = (maxlength - strlen(message)) / 2;
    printf("%*s%s%*s\n", padlen-1, "", message, padlen, "");
    partition(1, maxlength);
}

void option_header(char *names, size_t maxlength) {
    int padlen = (maxlength - strlen(names)) / 2;
    //partition(1, maxlength);
    printf("%*s%s%*s\n", padlen, "", names, padlen, "");
    partition(1, maxlength);
}

unsigned int choose_option(char **options, size_t n, size_t maxlength, char *category) {
    unsigned int choice;
    for (int i = 0; i < n; i++) {
        printf("%s%-3d%s%s\n", "  ", i+1, "- ", options[i]);
    }
    printf("%s%-3d%s\n", "  ", 0, "- EXIT");
    partition(1, maxlength);
    printf("%s%s%s","  Choose a ", category, " for the Analysis: ");
    scanf(" %u",&choice);
    
    return choice;
}

void mixed_header(char *names1, char *names2, size_t maxlength) {
    int padlen = (maxlength - strlen(names1) - strlen(names2)) / 2;
    printf("%*s%s%s%s%*s\n", padlen, "", names1, " - ", names2 , padlen, "");
    partition(1, maxlength);
}

void invalid_option(unsigned int option, char* category, size_t maxlength) {
    printf("%s%s%s%u%s", "  ", category, " ", option, " : Invalid option. Press any key to exit program...\n");
    partition(1, maxlength);
    while(getchar()!='\n'){}
    getchar(); // wait for ENTER
    exit(0);
}

void end_of_execution(size_t maxlength) {
    printf("%s", "\n\n  Execution ended successfully! Press any key to exit program...\n");
    partition(1, maxlength);
    //while(getchar()!='\n'){}
    getchar(); // wait for ENTER
}

void identify_simulation(unsigned int toolbox, unsigned int *system, unsigned int *module, char **toolboxesNames, char **systemNames, char** moduleNames, size_t numofsystems, size_t maxlength, size_t numofmodules) {
    clear_screen();
    welcome_header(maxlength);
    option_header(toolboxesNames[toolbox-1], maxlength);
    (*system) = choose_option(systemNames, numofsystems, maxlength, "System");
    // Chose System
    if (((*system) > 0) && ((*system) < numofsystems + 1)) {
        // Chose Module
        clear_screen();
        welcome_header(maxlength);
        option_header(toolboxesNames[toolbox - 1], maxlength);
        option_header(systemNames[(*system) - 1], maxlength);
        (*module) = choose_option(moduleNames, numofmodules, maxlength, "Module");
        if (((*module) > 0) && ((*module) < numofmodules + 1)) {
            clear_screen();
            welcome_header(maxlength);
            option_header(toolboxesNames[toolbox - 1], maxlength);
            mixed_header(systemNames[(*system) - 1], moduleNames[(*module) - 1], maxlength);
        }
        else if ((*module) == 0) {
            clear_screen();
            exit(0);
        }
        else {
            clear_screen();
            welcome_header(maxlength);
            option_header(toolboxesNames[toolbox - 1], maxlength);
            option_header(systemNames[(*system) - 1], maxlength);
            invalid_option((*module), "Module", maxlength);
        }
    }
    else if ((*system) == 0) {
        clear_screen();
        exit(0);
    }
    else {
        clear_screen();
        welcome_header(maxlength);
        option_header(toolboxesNames[toolbox - 1], maxlength);
        invalid_option((*system), "System", maxlength);
    }
}

int int_length(int value) {
    if (value == 0) {
        return 1;
    } 
    else {
        return (floor(log10(abs(value))) + 1);
    }
}

void progress_bar_old(int mode, double var, double var_i, double var_f) {
    double perc;
    // Actual percentage
    if (mode == 1) {
        perc = 100;
    }
    else {
        perc = (var/(var_f - var_i))*100;
    }
    // Filled Part of the progress bar
    int fill = (perc * PROGRESS_BAR_LEN) / 100;  
    printf("\r  Progress: |");
    for(int i = 0; i < fill; i++) {
        printf("#");
    }
    // Unfilled part of the progress bar
    for (int i = 0; i < PROGRESS_BAR_LEN - fill; i++) {
        printf(".");
    }
    if (perc > 100) {
        perc = 100;
    }
    printf("| %.1lf %% ", perc);
    fflush(stdout);
}

static void print_progress_bar(double perc) {
    // Filled Part of the progress bar
    int fill = (perc * PROGRESS_BAR_LEN) / 100;  
    printf("\r  Progress: |");
    for(int i = 0; i < fill; i++) {
        printf("#");
    }
    // Unfilled part of the progress bar
    for (int i = 0; i < PROGRESS_BAR_LEN - fill; i++) {
        printf(".");
    }
    if (perc > 100) {
        perc = 100;
    }
    printf("| %.1lf %% ", perc);
    fflush(stdout);
}

void progress_bar(int mode, double var, double var_i, double var_f) {
    double perc;
    // Actual percentage
    if (mode == 1) {
        perc = 100;
    }
    else {
        perc = (var/(var_f - var_i))*100;
    }
    // Print bar
    if (perc < 33) {
        red();
        print_progress_bar(perc);
        reset_color();
    }
    else if ((perc >= 33) && (perc < 66)) {
        yellow();
        print_progress_bar(perc);
        reset_color();
    }
    else if (perc >= 66) {
        green();
        print_progress_bar(perc);
        reset_color();
    }
    else {
        print_progress_bar(perc);
    }
}

// Simulation Prints

void fpartition(FILE *output_file, int mode, size_t maxlength) {
    if (mode == 1) {
        fprintf(output_file, " ");
        for (size_t i = 0; i < maxlength-2; i++) {
            fprintf(output_file, "=");
        }
        fprintf(output_file, " ");
    }
    else if (mode == 2) {
        fprintf(output_file, " ");
        for (size_t i = 0; i < maxlength-2; i++) {
            fprintf(output_file, "-");
        }
        fprintf(output_file, " ");
    }
    else {
        printf("DEBUG WARNING: Please select 1 for main partition or 2 for secondary partition for printing info in the output file");
        exit(0);
    }
    fprintf(output_file, "\n");
}

void print_list_of_indexes(int n, int *indexes, int spcvalue, int spcname, char* name) {
    // Allocate memory to handle operations    
    char *buffer = calloc((spcvalue + 1), sizeof(char));  // Buffer to store each index at a time
    char *newrow = calloc((spcvalue + 1), sizeof(char)); // String to store each row
    char *formatedrow = calloc((spcvalue + spcname + 2), sizeof(char)); // String to store each formatted row
    
    // Declare counter for number of rows found
    int rows = 0; 
    // Sweep across all indexes
    for (int i = 0; i < n; i++) {
        // Check if is the last index
        if (i == n - 1) {
            // Put the final index into the buffer
            sprintf(buffer, "%d", indexes[i]);
            // Check if the length of buffer and newrow together are bigger than the maximum space for values
            if (strlen(newrow) + strlen(buffer) > spcvalue) {
                // Format the row with names
                if (rows == 0) {
                    sprintf(formatedrow, "%-*s %-*s", spcname, name, spcvalue, newrow);
                }
                else {
                    sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
                }
                printf("\n%s", formatedrow);
                // Empty the row string
                strcpy(newrow, "\0");
                strcpy(formatedrow, "\0");
            }
            // Update newrow with the new value
            strcat(newrow, buffer);
            // Format the row with names
            if (rows == 0) {
                sprintf(formatedrow, "%-*s %-*s", spcname, name, spcvalue, newrow);
                printf("%s", formatedrow);    
            }
            else {
                sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
                printf("\n%s", formatedrow);
            }
        }
        else {
            // Put the next index into the buffer
            sprintf(buffer, "%d, ", indexes[i]);
            // Check if the length of buffer and newrow together are bigger than the maximum space for values
            if (strlen(newrow) + strlen(buffer) > spcvalue) {
                rows = rows + 1;
                // Format the row with names
                if (rows == 1) {
                    sprintf(formatedrow, "%-*s %-*s", spcname, name, spcvalue, newrow);
                    // Print formatted row
                    printf("%s", formatedrow);
                }
                else {
                    sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
                    // Print formatted row
                    printf("\n%s", formatedrow);
                }
                // Empty the row string
                strcpy(newrow, "\0");
                strcpy(formatedrow, "\0");
            }
            // Update newrow with the new value
            strcat(newrow, buffer);
        }
    }
    // Free Memory
    free(newrow); free(buffer); free(formatedrow);
}

void fprint_list_of_indexes(FILE *output_file, int n, int *indexes, int spcvalue, int spcname, char* name) {
    // Allocate memory to handle operations
    char *buffer = calloc((spcvalue + 1), sizeof(char));  // Buffer to store each index at a time
    char *newrow = calloc((spcvalue + 1), sizeof(char)); // String to store each row
    char *formatedrow = calloc((spcvalue + spcname + 2), sizeof(char)); // String to store each formatted row
    // Declare counter for number of rows found
    int rows = 0; 
    // Sweep across all indexes
    for (int i = 0; i < n; i++) {
        // Check if is the last index
        if (i == n - 1) {
            // Put the final index into the buffer
            sprintf(buffer, "%d", indexes[i]);
            // Check if the length of buffer and newrow together are bigger than the maximum space for values
            if (strlen(newrow) + strlen(buffer) > spcvalue) {
                // Format the row with names
                if (rows == 0) {
                    sprintf(formatedrow, "%-*s %-*s", spcname, name, spcvalue, newrow);
                }
                else {
                    sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
                }
                fprintf(output_file, "\n%s", formatedrow);
                // Empty the row string
                strcpy(newrow, "");
                strcpy(formatedrow, "");
            }
            // Update newrow with the new value
            strcat(newrow, buffer);
            // Format the row with names
            if (rows == 0) {
                sprintf(formatedrow, "%-*s %-*s", spcname, name, spcvalue, newrow);
                fprintf(output_file, "%s", formatedrow);    
            }
            else {
                sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
                fprintf(output_file, "\n%s", formatedrow);
            }
        }
        else {
            // Put the next index into the buffer
            sprintf(buffer, "%d, ", indexes[i]);
            // Check if the length of buffer and newrow together are bigger than the maximum space for values
            if (strlen(newrow) + strlen(buffer) > spcvalue) {
                rows = rows + 1;
                // Format the row with names
                if (rows == 1) {
                    sprintf(formatedrow, "%-*s %-*s", spcname, name, spcvalue, newrow);
                    // Print formatted row
                    fprintf(output_file, "%s", formatedrow);
                }
                else {
                    sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
                    // Print formatted row
                    fprintf(output_file, "\n%s", formatedrow);
                }
                // Empty the row string
                strcpy(newrow, "");
                strcpy(formatedrow, "");
            }
            // Update newrow with the new value
            strcat(newrow, buffer);
        }
    }
    // Free Memory
    free(newrow); free(buffer); free(formatedrow);
}

void write_initial_conditions(int dim, double *x, double t, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  Initial Conditions\n");
    partition(2, maxlength);
    printf("%-*s %-*.15lf\n", spcname, "  Initial Time (t):", spcvalue, t);
    for (int i = 0; i < dim; i++) {
        printf("%s%d%-*s %-*.15lf\n", "  x[", i, spcname - 4 - int_length(i),"]:", spcvalue, x[i]);
    }
}

void fwrite_initial_conditions(FILE *output_file, int dim, double *x, double t, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  Initial Conditions\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*.15lf\n", spcname, "  Initial Time (t):", spcvalue, t);
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "%s%d%-*s %-*.15lf\n", "  x[", i, spcname - 4 - int_length(i),"]:", spcvalue, x[i]);
    }
}

void write_sys_parameters(int npar, double *par, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  System Parameters\n");
    partition(2, maxlength);
    for (int i = 0; i < npar; i++) {
        printf("%s%d%-*s %-*g\n", "  par[", i, spcname - 6 - int_length(i), "]:", spcvalue, par[i]);
    }
}

void fwrite_sys_parameters(FILE *output_file, int npar, double *par, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  System Parameters\n");
    fpartition(output_file, 2, maxlength);
    for (int i = 0; i < npar; i++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  par[", i, spcname - 6 - int_length(i), "]:", spcvalue, par[i]);
    }
}

void write_RMS_calculations_info(int n, int *index, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  RMS Calculation Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Number of RMS Calculations:", spcvalue, n);
    print_list_of_indexes(n, index, spcvalue, spcname, "  State Variables Indexes:");
    if (n > 0) {
        printf("\n");
    }
}

void fwrite_RMS_calculations_info(FILE *output_file, int n, int *index, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  RMS Calculation Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of RMS Calculations:", spcvalue, n);
    fprint_list_of_indexes(output_file, n, index, spcvalue, spcname, "  State Variables Indexes:");
    if (n > 0) {
        fprintf(output_file, "\n");
    }
}

void write_custom_info_calculations(int n, int nf, int *findex, int nscr, int *scrindex, size_t maxlength, double percname) {
    if (n > 0) {
        int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
        int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
        partition(2, maxlength);
        printf("  Custom Calculation Parameters\n");
        partition(2, maxlength);
        printf("%-*s %-*d\n", spcname, "  Number of Custom Calculations:", spcvalue, n);
        if (nf > 0) {
            printf("%-*s %-*d\n", spcname, "  Custom Values Printed on File:", spcvalue, nf);
            
            print_list_of_indexes(nf, findex, spcvalue, spcname, "  Printed on File Indexes:");
        }
        if (nscr > 0) {
            printf("\n");
            printf("%-*s %-*d\n", spcname, "  Custom Values Printed on Screen:", spcvalue, nscr);
            print_list_of_indexes(nscr, scrindex, spcvalue, spcname, "  Printed on Screen Indexes:");
        }
        printf("\n");
    }
}

void fwrite_custom_info_calculations(FILE* output_file, int n, int nf, int *findex, int nscr, int *scrindex, size_t maxlength, double percname) {
    if (n > 0) {
        int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
        int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
        fpartition(output_file, 2, maxlength);
        fprintf(output_file, "  Custom Calculation Parameters\n");
        fpartition(output_file, 2, maxlength);
        fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Custom Calculations:", spcvalue, n);
        if (nf > 0) {
            fprintf(output_file, "%-*s %-*d\n", spcname, "  Custom Values Printed on File:", spcvalue, nf);
            fprint_list_of_indexes(output_file, nf, findex, spcvalue, spcname, "  Printed on File Indexes:");
        }
        if (nscr > 0) {
            fprintf(output_file, "\n");
            fprintf(output_file, "%-*s %-*d\n", spcname, "  Custom Values Printed on Screen:", spcvalue, nscr);
            fprint_list_of_indexes(output_file, nscr, scrindex, spcvalue, spcname, "  Printed on Screen Indexes:");
        }
        fprintf(output_file, "\n");
    }
}

void fwrite_module_and_system(FILE *output_file, char *funcname, char *modulename, size_t maxlength) {
    time_t tm;
    time(&tm);
    fprintf(output_file, "  Date/Time:  %s", ctime(&tm)); 
    fpartition(output_file, 1, maxlength);
    fprintf(output_file, "  %s: %s\n", modulename, funcname);
    fpartition(output_file, 1, maxlength);;
}

void write_prog_parameters_timeseries(int dim, int npar, int np, int ndiv, int trans, double h, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*g\n", spcname, "  Timestep Value:", spcvalue, h);
}

void fwrite_prog_parameters_timeseries(FILE* output_file, char *funcname, int dim, int npar, int np, int ndiv, int trans, double h, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    //Get time and date
    fwrite_module_and_system(output_file, funcname, "Time Series", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Timestep Value:", spcvalue, h);
}

void write_prog_parameters_ftimeseries(int dim, int npar, int maxper, int np, int ndiv, int trans, double h, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*g\n", spcname, "  Timestep Value:", spcvalue, h);
    printf("%-*s %-*d\n", spcname, "  Max Periodicity Considered:", spcvalue, maxper);
}

void fwrite_prog_parameters_ftimeseries(FILE* output_file, char *funcname, int dim, int npar, int maxper, int np, int ndiv, int trans, double h, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Time Series", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Timestep Value:", spcvalue, h);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Max Periodicity Considered:", spcvalue, maxper);
}

void write_prog_parameters_strobopoinc(int dim, int npar, int np, int ndiv, int trans, double h, size_t maxlength,  double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*g\n", spcname, "  Timestep Value:", spcvalue, h);
}

void fwrite_prog_parameters_strobopoinc(FILE* output_file, char *funcname, int dim, int npar, int np, int ndiv, int trans, double h, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Stroboscopic Poincaré Map", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Timestep Value:", spcvalue, h);
}

void write_prog_parameters_lyapunov(int dim, int npar, int np, int ndiv, int trans, double h, size_t maxlength,  double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*g\n", spcname, "  Timestep Value:", spcvalue, h);
}

void fwrite_prog_parameters_lyapunov(FILE *output_file, char *funcname, int dim, int npar, int np, int ndiv, int trans, double h, size_t maxlength,  double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Lyapunov Exponents (Wolf)", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Timestep Value:", spcvalue, h);
}

void write_prog_parameters_bifurcation(int dim, int npar, int np, int ndiv, int trans, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void fwrite_prog_parameters_bifurcation(FILE* output_file, char *funcname, int dim, int npar, int np, int ndiv, int trans, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength;  // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;       // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Bifurcation Diagram", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void write_bifurcation_info(double *parrange, int parindex, int bifmode, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength;  // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;       // Space for value of variable
    partition(2, maxlength);
    printf("  Bifurcation Parameters\n");
    partition(2, maxlength);
    if (bifmode == 0) {
        printf("%-*s %-*s\n", spcname, "  Bifucation Mode:", spcvalue, "Following Attractor");
    }
    else if (bifmode == 1) {
        printf("%-*s %-*s\n", spcname, "  Bifucation Mode:", spcvalue, "Reseting Initial Conditions");
    } 
    else {
        printf("\n  Invalid Bifurcation Mode of %d...\nChange to 0 or 1 in the input file and run the program again.\nExiting Program...!\n", bifmode);
        exit(1);
    }
    printf("%-*s %-*d\n", spcname, "  Control Parameter Index:", spcvalue, parindex);
    printf("%-*s %-*g\n", spcname, "  Initial Control Parameter:", spcvalue, parrange[0]);
    printf("%-*s %-*g\n", spcname, "  Final Control Parameter:", spcvalue, parrange[1]);
    printf("%-*s %-*g\n", spcname, "  Control Parameter Step:", spcvalue, (parrange[1] - parrange[0]) / (parrange[2] - 1));
    printf("%-*s %-*g\n", spcname, "  Total Number of Steps:", spcvalue, parrange[2]);
}

void fwrite_bifurcation_info(FILE* output_file, double *parrange, int parindex, int bifmode, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength;  // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;       // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  Bifurcation Parameters\n");
    fpartition(output_file, 2, maxlength);
    
    if (bifmode == 0) {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Bifucation Mode:", spcvalue, "Following Attractor");
    }
    else if (bifmode == 1) {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Bifucation Mode:", spcvalue, "Reseting Initial Conditions");
    } 
    else {
        fprintf(output_file, "\n  Invalid Bifurcation Mode of %d...\nChange to 0 or 1 in the input file and run the program again.\nExiting Program...!\n", bifmode);
        exit(1);
    }
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Control Parameter Index:", spcvalue, parindex);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Initial Control Parameter:", spcvalue, parrange[0]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Final Control Parameter:", spcvalue, parrange[1]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Control Parameter Step:", spcvalue, (parrange[1] - parrange[0]) / (parrange[2] - 1));
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Total Number of Steps:", spcvalue, parrange[2]);
}

void write_prog_parameters_fbifurcation(int dim, int npar, int np, int ndiv, int trans, int maxper, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*d\n", spcname, "  Max Periodicity Considered:", spcvalue, maxper);
    printf("%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void fwrite_prog_parameters_fbifurcation(FILE *output_file, char *funcname, int dim, int npar, int np, int ndiv, int trans, int maxper, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Bifurcation Diagram", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Max Periodicity Considered:", spcvalue, maxper);
    fprintf(output_file, "%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void write_prog_parameters_dyndiag(int dim, int npar, int np, int ndiv, int trans, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void fwrite_prog_parameters_dyndiag(FILE *output_file, char* funcname, int dim, int npar, int np, int ndiv, int trans, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Dynamical Response Diagram", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void write_dyndiag_info(double *parrange, int indexX, int indexY, int bifmode, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength;  // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;       // Space for value of variable
    partition(2, maxlength);
    printf("  Dynamical Diagram Parameters\n");
    partition(2, maxlength);
    if (bifmode == 0) {
        printf("%-*s %-*s\n", spcname, "  Diagram Mode:", spcvalue, "Following Attractor");
    }
    else if (bifmode == 1) {
        printf("%-*s %-*s\n", spcname, "  Diagram Mode:", spcvalue, "Reseting Initial Conditions");
    } 
    else {
        printf("\n  Invalid Diagram Mode of %d...\nChange to 0 or 1 in the input file and run the program again.\nExiting Program...!\n", bifmode);
        exit(1);
    }
    printf("%-*s %g x %-*g\n", spcname, "  Resolution:", parrange[2], spcvalue - 3 - int_length((int)parrange[2]), parrange[5]);
    printf("\n");
    printf("%-*s %-*d\n", spcname, "  Control Parameter Index (x):", spcvalue, indexX);
    printf("%-*s %-*g\n", spcname, "  Initial Control Parameter (x):", spcvalue, parrange[0]);
    printf("%-*s %-*g\n", spcname, "  Final Control Parameter (x):", spcvalue, parrange[1]);
    printf("%-*s %-*g\n", spcname, "  Control Parameter Step (x):", spcvalue, (parrange[1] - parrange[0]) / (parrange[2] - 1));
    printf("%-*s %-*g\n", spcname, "  Total Number of Steps (x):", spcvalue, parrange[2]);
    printf("\n");
    printf("%-*s %-*d\n", spcname, "  Control Parameter Index (y):", spcvalue, indexY);
    printf("%-*s %-*g\n", spcname, "  Initial Control Parameter (y):", spcvalue, parrange[3]);
    printf("%-*s %-*g\n", spcname, "  Final Control Parameter (y):", spcvalue, parrange[4]);
    printf("%-*s %-*g\n", spcname, "  Control Parameter Step (y):", spcvalue, (parrange[4] - parrange[3]) / (parrange[5] - 1));
    printf("%-*s %-*g\n", spcname, "  Total Number of Steps (y):", spcvalue, parrange[5]);
}

void fwrite_dyndiag_info(FILE *output_file, double *parrange, int indexX, int indexY, int bifmode, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength;  // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;       // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  Dynamical Diagram Parameters\n");
    fpartition(output_file, 2, maxlength);
    if (bifmode == 0) {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Diagram Mode:", spcvalue, "Following Attractor");
    }
    else if (bifmode == 1) {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Diagram Mode:", spcvalue, "Reseting Initial Conditions");
    } 
    else {
        fprintf(output_file, "\n  Invalid Diagram Mode of %d...\nChange to 0 or 1 in the input file and run the program again.\nExiting Program...!\n", bifmode);
        exit(1);
    }
    fprintf(output_file, "%-*s %g x %-*g\n", spcname, "  Resolution:", parrange[2], spcvalue - 3 - int_length((int)parrange[2]), parrange[5]);
    fprintf(output_file, "\n");
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Control Parameter Index (x):", spcvalue, indexX);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Initial Control Parameter (x):", spcvalue, parrange[0]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Final Control Parameter (x):", spcvalue, parrange[1]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Control Parameter Step (x):", spcvalue, (parrange[1] - parrange[0]) / (parrange[2] - 1));
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Total Number of Steps (x):", spcvalue, parrange[2]);
    fprintf(output_file, "\n");
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Control Parameter Index (y):", spcvalue, indexY);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Initial Control Parameter (y):", spcvalue, parrange[3]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Final Control Parameter (y):", spcvalue, parrange[4]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Control Parameter Step (y):", spcvalue, (parrange[4] - parrange[3]) / (parrange[5] - 1));
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Total Number of Steps (y):", spcvalue, parrange[5]);
}

void write_prog_parameters_fdyndiag(int dim, int npar, int np, int ndiv, int maxper, int trans, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*d\n", spcname, "  Max Periodicity Considered:", spcvalue, maxper);
    printf("%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void fwrite_prog_parameters_fdyndiag(FILE *output_file, char *funcname, int dim, int npar, int np, int ndiv, int maxper, int trans, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Full Dynamical Response Diagram", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Max Periodicity Considered:", spcvalue, maxper);
    fprintf(output_file, "%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void write_prog_parameters_fforcbasin(int dim, int npar, int np, int ndiv, int maxper, int trans, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    printf("%-*s %-*d\n", spcname, "  Max Periodicity Considered:", spcvalue, maxper);
    printf("%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void fwrite_prog_parameters_fforcbasin(FILE *output_file, char *funcname, int dim, int npar, int np, int ndiv, int maxper, int trans, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Full Basin of Attraction (Forced)", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Transient Considered:", spcvalue, trans);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Max Periodicity Considered:", spcvalue, maxper);
    fprintf(output_file, "%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void write_prog_parameters_epbasin(int dim, int npar, int np, int ndiv, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("\n  Program Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    printf("%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    printf("%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    printf("%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    printf("%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    printf("%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void fwrite_prog_parameters_epbasin(FILE *output_file, char *funcname, int dim, int npar, int np, int ndiv, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fwrite_module_and_system(output_file, funcname, "Basin of Attraction (Fixed Points)", maxlength);
    fprintf(output_file, "\n  Program Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Dimension:", spcvalue, dim);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of Parameters:", spcvalue, npar);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Forcing Periods:", spcvalue, np);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Timesteps per Period:", spcvalue, ndiv);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Total Number of Timesteps:", spcvalue, np*ndiv);
    fprintf(output_file, "%-*s %-*s\n", spcname, "  Timestep Value:", spcvalue, "(2*pi)/(nDiv*par[0])");
}

void write_basin_info(char *type, double *icrange, int indexX, int indexY, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength;  // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;       // Space for value of variable
    partition(2, maxlength);
    printf("  Basin of Attraction (%s) Parameters\n", type);
    partition(2, maxlength);
    printf("%-*s %g x %-*g\n", spcname, "  Resolution:", icrange[2], spcvalue - 3 - int_length((int)icrange[2]), icrange[5]);
    printf("\n");
    printf("%-*s %-*d\n", spcname, "  Initial Condition Index (x):", spcvalue, indexX);
    printf("%-*s %-*g\n", spcname, "  Initial Value (x):", spcvalue, icrange[0]);
    printf("%-*s %-*g\n", spcname, "  Final Value (x):", spcvalue, icrange[1]);
    printf("%-*s %-*g\n", spcname, "  Initial Contition Step (x):", spcvalue, (icrange[1] - icrange[0]) / (icrange[2] - 1));
    printf("%-*s %-*g\n", spcname, "  Total Number of Steps (x):", spcvalue, icrange[2]);
    printf("\n");
    printf("%-*s %-*d\n", spcname, "  Initial Condition Index (y):", spcvalue, indexY);
    printf("%-*s %-*g\n", spcname, "  Initial Value (y):", spcvalue, icrange[3]);
    printf("%-*s %-*g\n", spcname, "  Final Value (y):", spcvalue, icrange[4]);
    printf("%-*s %-*g\n", spcname, "  Initial Condition Step (y):", spcvalue, (icrange[4] - icrange[3]) / (icrange[5] - 1));
    printf("%-*s %-*g\n", spcname, "  Total Number of Steps (y):", spcvalue, icrange[5]);
}

void fwrite_basin_info(char* type, FILE *output_file, double *icrange, int indexX, int indexY, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength;  // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;       // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  Basin of Attraction (%s) Parameters\n", type);
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %g x %-*g\n", spcname, "  Resolution:", icrange[2], spcvalue - 3 - int_length((int)icrange[2]), icrange[5]);
    fprintf(output_file, "\n");
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Initial Condition Index (x):", spcvalue, indexX);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Initial Value (x):", spcvalue, icrange[0]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Final Value (x):", spcvalue, icrange[1]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Initial Contition Step (x):", spcvalue, (icrange[1] - icrange[0]) / (icrange[2] - 1));
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Total Number of Steps (x):", spcvalue, icrange[2]);
    fprintf(output_file, "\n");
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Initial Condition Index (y):", spcvalue, indexY);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Initial Value (y):", spcvalue, icrange[3]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Final Value (y):", spcvalue, icrange[4]);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Initial Condition Step (y):", spcvalue, (icrange[4] - icrange[3]) / (icrange[5] - 1));
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Total Number of Steps (y):", spcvalue, icrange[5]);
}

void print_attractor(int attrac, int maxper, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    // Get the attractor and convert into a string
    char number[5];
    sprintf(number, "%d", attrac);

    if (attrac < maxper) { 
        printf("%-*s %s%d%-*s\n", spcname, "  Type of Motion:", "Period-", attrac, (int)(spcvalue - 8 - strlen(number) - 1), "T");
    }
    else if (attrac == maxper) {
        printf("%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Many Periods");
    }
    else if (attrac == maxper + 1) {
        printf("%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Chaotic");
    }
    else if (attrac == maxper + 2) {
        printf("%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Hyperchaotic");
    }
    else if (attrac == maxper + 3) {
        printf("%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Undefined");
    }
    else {
        printf("%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Undefined (Escape)");
    } 
    partition(2, maxlength);
}

void fprint_attractor(FILE *output_file, int attrac, int maxper, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    // Get the attractor and convert into a string
    char number[5];
    sprintf(number, "%d", attrac);

    if (attrac < maxper) { 
        fprintf(output_file, "%-*s %s%d%-*s\n", spcname, "  Type of Motion:", "Period-", attrac, (int)(spcvalue - 8 - strlen(number) - 1), "T");
    }
    else if (attrac == maxper) {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Many Periods");
    }
    else if (attrac == maxper + 1) {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Chaotic");
    }
    else if (attrac == maxper + 2) {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Hyperchaotic");
    }
    else if (attrac == maxper + 3) {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Undefined");
    }
    else {
        fprintf(output_file, "%-*s %-*s\n", spcname, "  Type of Motion:", spcvalue, "Undefined (Escape)");
    } 
    fpartition(output_file, 2, maxlength);
}

void print_RMS(int nRMS, int *rmsindex, double *xRMS, double *overallxRMS, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  Results of RMS Calculations\n");
    partition(2, maxlength);
    // Print RMS calculations on screen
    for (int q = 0; q < nRMS; q++) {
        printf("%s%d%-*s %-*g\n", "  xRMS[", rmsindex[q], spcname - 7 - int_length(rmsindex[q]), "]:", spcvalue, xRMS[rmsindex[q]]);
    }
    for (int q = 0; q < nRMS; q++) {
        printf("%s%d%-*s %-*g\n", "  Overall xRMS[", rmsindex[q], spcname - 15 - int_length(rmsindex[q]), "]:", spcvalue, overallxRMS[rmsindex[q]]);
    }
}

void fprint_RMS(FILE *output_file, int nRMS, int *rmsindex, double *xRMS, double *overallxRMS, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  Results of RMS Calculations\n");
    fpartition(output_file, 2, maxlength);
    // Print RMS calculations on file
    for (int q = 0; q < nRMS; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  xRMS[", rmsindex[q], spcname - 7 - int_length(rmsindex[q]), "]:", spcvalue, xRMS[rmsindex[q]]);
    }
    for (int q = 0; q < nRMS; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  Overall xRMS[", rmsindex[q], spcname - 15 - int_length(rmsindex[q]), "]:", spcvalue, overallxRMS[rmsindex[q]]);
    }
}

void print_customcalc(int nprintscr, int *printscrindex, double *customvalue, char **customnames, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  Results of Custom Calculations\n");
    partition(2, maxlength);
    // Print Custom calculations on screen
    for (int q = 0; q < nprintscr; q++) {
        printf("  %s%-*s %-*g\n", customnames[printscrindex[q]], (int)(spcname - strlen(customnames[printscrindex[q]]) - 2), ":", spcvalue, customvalue[printscrindex[q]]);
    }
}

void fprint_customcalc(FILE *output_file, int nprintscr, int *printscrindex, double *customvalue, char **customnames, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  Results of Custom Calculations\n");
    fpartition(output_file, 2, maxlength);
    // Print Custom calculations on screen
    for (int q = 0; q < nprintscr; q++) {
        fprintf(output_file, "  %s%-*s %-*g\n", customnames[printscrindex[q]], (int)(spcname - strlen(customnames[printscrindex[q]]) - 2), ":", spcvalue, customvalue[printscrindex[q]]);
    }
}

void print_minmax(double *xmin, double *xmax, double *overallxmin, double *overallxmax, int dim, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  Minimum and Maximum Values\n");
    partition(2, maxlength);
    // Print min and maximum values
    for (int q = 0; q < dim; q++) {
        printf("%s%d%-*s %-*g\n", "  xmin[", q, spcname - 7 - int_length(q), "]:", spcvalue, xmin[q]);
    }
    for (int q = 0; q < dim; q++) {
        printf("%s%d%-*s %-*g\n", "  xmax[", q, spcname - 7 - int_length(q), "]:", spcvalue, xmax[q]);
    }
    for (int q = 0; q < dim; q++) {
        printf("%s%d%-*s %-*g\n", "  Overall xmin[", q, spcname - 15 - int_length(q), "]:", spcvalue, overallxmin[q]);
    }
    for (int q = 0; q < dim; q++) {
        printf("%s%d%-*s %-*g\n", "  Overall xmax[", q, spcname - 15 - int_length(q), "]:", spcvalue, overallxmax[q]);
    }
}

void fprint_minmax(FILE *output_file, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int dim, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  Minimum and Maximum Values\n");
    fpartition(output_file, 2, maxlength);
    // Print min and maximum values
    for (int q = 0; q < dim; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  xMIN[", q, spcname - 7 - int_length(q), "]:", spcvalue, xmin[q]);
    }
    for (int q = 0; q < dim; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  xMAX[", q, spcname - 7 - int_length(q), "]:", spcvalue, xmax[q]);
    }
    for (int q = 0; q < dim; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  Overall xMIN[", q, spcname - 15 - int_length(q), "]:", spcvalue, overallxmin[q]);
    }
    for (int q = 0; q < dim; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  Overall xMAX[", q, spcname - 15 - int_length(q), "]:", spcvalue, overallxmax[q]);
    }
}

