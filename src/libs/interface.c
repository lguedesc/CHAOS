#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void clear_screen() {
    #ifdef _WIN32
        system("cls");
    #else
        system("clear");
    #endif
}

// Main Interface 
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
    char *message = "Welcome to CHAOS, a Nonlinear Dynamics Package for Harmonically Forced Systems";
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
    while(getchar()!='\n'){}
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

void write_initial_conditions(int dim, double *x, double t, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  Initial Conditions\n");
    partition(2, maxlength);
    printf("%-*s %-*g\n", spcname, "  Initial Time (t):", spcvalue, t);
    for (int i = 0; i < dim; i++) {
        printf("%s%d%-*s %-*g\n", "  x[", i, spcname - 5,"]:", spcvalue, x[i]);
    }
}

void fwrite_initial_conditions(FILE *output_file, int dim, double *x, double t, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  Initial Conditions\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*g\n", spcname, "  Initial Time (t):", spcvalue, t);
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  x[", i, spcname - 5,"]:", spcvalue, x[i]);
    }
}

void write_sys_parameters(int npar, double *par, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  System Parameters\n");
    partition(2, maxlength);
    for (int i = 0; i < npar; i++) {
        printf("%s%d%-*s %-*g\n", "  par[", i, spcname - 7, "]:", spcvalue, par[i]);
    }
}

void fwrite_sys_parameters(FILE *output_file, int npar, double *par, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  System Parameters\n");
    fpartition(output_file, 2, maxlength);
    for (int i = 0; i < npar; i++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  par[", i, spcname - 7, "]:", spcvalue, par[i]);
    }
}

void write_RMS_calculations_info(int n, int *index, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  RMS Calculation Parameters\n");
    partition(2, maxlength);
    printf("%-*s %-*d\n", spcname, "  Number of RMS Calculations:", spcvalue, n);
    printf("%-*s %d", spcname, "  State Variables Indexes:", index[0]);
    for (int i = 1; i < n; i++) {
        printf(", %d", index[i]);
    }
    printf("\n");
}

void fwrite_RMS_calculations_info(FILE *output_file, int n, int *index, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "  RMS Calculation Parameters\n");
    fpartition(output_file, 2, maxlength);
    fprintf(output_file, "%-*s %-*d\n", spcname, "  Number of RMS Calculations:", spcvalue, n);
    fprintf(output_file, "%-*s %d", spcname, "  State Variables Indexes:", index[0]);
    for (int i = 1; i < n; i++) {
        fprintf(output_file, ", %d", index[i]);
    }
    fprintf(output_file, "\n");
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
            printf("%-*s %d", spcname, "  Printed on File Indexes:", findex[0]);
            for (int i = 1; i < nf; i++) {
                printf(", %d", findex[i]);
            }
        }
        if (nscr > 0) {
            printf("\n");
            printf("%-*s %-*d\n", spcname, "  Custom Values Printed on Screen:", spcvalue, nscr);
            printf("%-*s %d", spcname, "  Printed on Screen Indexes:", scrindex[0]);
            for (int i = 1; i < nscr; i++) {
                printf(", %d", scrindex[i]);
            }
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
            fprintf(output_file, "%-*s %d", spcname, "  Printed on File Indexes:", findex[0]);
            for (int i = 1; i < nf; i++) {
                fprintf(output_file, ", %d", findex[i]);
            }
        }
        if (nscr > 0) {
            fprintf(output_file, "\n");
            fprintf(output_file, "%-*s %-*d\n", spcname, "  Custom Values Printed on Screen:", spcvalue, nscr);
            fprintf(output_file, "%-*s %d", spcname, "  Printed on Screen Indexes:", scrindex[0]);
            for (int i = 1; i < nscr; i++) {
                fprintf(output_file, ", %d", scrindex[i]);
            }
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

void print_attractor(int attrac, int maxper, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    printf("  Attractor\n");
    partition(2, maxlength);
    // Get the attractor and convert into a string
    char number[5];
    sprintf(number, "%d", attrac);

    if (attrac < maxper) { 
        printf("%-*s %s%d%-*s\n", spcname, "  Type of Motion:", "Period-", attrac, spcvalue - 8 - strlen(number) - 1, "T");
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
}

void fprint_attractor(FILE *output_file, int attrac, int maxper, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    fprintf(output_file, "  Attractor\n");
    fpartition(output_file, 2, maxlength);
    // Get the attractor and convert into a string
    char number[5];
    sprintf(number, "%d", attrac);

    if (attrac < maxper) { 
        fprintf(output_file, "%-*s %s%d%-*s\n", spcname, "  Type of Motion:", "Period-", attrac, spcvalue - 8 - strlen(number) - 1, "T");
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
}

void print_RMS(int nRMS, int *rmsindex, double *xRMS, double *overallxRMS, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  Results of RMS Calculations\n");
    partition(2, maxlength);
    // Print RMS calculations on screen
    for (int q = 0; q < nRMS; q++) {
        printf("%s%d%-*s %-*g\n", "  xRMS[", rmsindex[q], spcname - 8, "]:", spcvalue, xRMS[rmsindex[q]]);
    }
    for (int q = 0; q < nRMS; q++) {
        printf("%s%d%-*s %-*g\n", "  Overall xRMS[", rmsindex[q], spcname - 16, "]:", spcvalue, overallxRMS[rmsindex[q]]);
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
        fprintf(output_file, "%s%d%-*s %-*g\n", "  xRMS[", rmsindex[q], spcname - 8, "]:", spcvalue, xRMS[rmsindex[q]]);
    }
    for (int q = 0; q < nRMS; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  Overall xRMS[", rmsindex[q], spcname - 16, "]:", spcvalue, overallxRMS[rmsindex[q]]);
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
        printf("%s%d%-*s %-*g\n", "  xmix[", q, spcname - 8, "]:", spcvalue, xmin[q]);
    }
    for (int q = 0; q < dim; q++) {
        printf("%s%d%-*s %-*g\n", "  xmax[", q, spcname - 8, "]:", spcvalue, xmax[q]);
    }
    for (int q = 0; q < dim; q++) {
        printf("%s%d%-*s %-*g\n", "  Overall xmin[", q, spcname - 16, "]:", spcvalue, overallxmin[q]);
    }
    for (int q = 0; q < dim; q++) {
        printf("%s%d%-*s %-*g\n", "  Overall xmax[", q, spcname - 16, "]:", spcvalue, overallxmax[q]);
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
        fprintf(output_file, "%s%d%-*s %-*g\n", "  xMIN[", q, spcname - 8, "]:", spcvalue, xmin[q]);
    }
    for (int q = 0; q < dim; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  xMAX[", q, spcname - 8, "]:", spcvalue, xmax[q]);
    }
    for (int q = 0; q < dim; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  Overall xMIN[", q, spcname - 16, "]:", spcvalue, overallxmin[q]);
    }
    for (int q = 0; q < dim; q++) {
        fprintf(output_file, "%s%d%-*s %-*g\n", "  Overall xMAX[", q, spcname - 16, "]:", spcvalue, overallxmax[q]);
    }
}

