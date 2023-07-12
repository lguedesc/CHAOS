#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <sys/stat.h>
#include "defines.h"
#include "msg.h"

#ifdef _WIN32
    #include <direct.h>
    const char SEP = '\\';
#else
    const char SEP = '/';
#endif    

// Define function to print name of the variable in a string
#define getName(var)  #var

struct stat info;
/*
char *convert_dir(const char *rawdir) {
    if (rawdir == NULL) {
        perror(rawdir);
        return NULL;
    }
    // Insert variables
    int n = 200;
    char *newdir = malloc(n*sizeof(char));
    char *ch;
    // Copy string pathname into string path
    strcpy(newdir, rawdir);
    // Read the string char by char and do the conversion based on OS
    for (ch = newdir+1; *ch; ch++) {
        #ifdef _WIN32
            if (*ch == '/') {
                *ch = '\\';
            }
        #else
            if (*ch == '\\') {
                *ch = '/';
            }
        #endif
    }
    // The user is responsible to free (newdir) after function call
    return newdir;
}
*/

char *convert_dir(const char *rawdir) {
    // Check for string
    if (rawdir == NULL) {
        perror("rawdir is NULL");
        return NULL;
    }
    // Get the length of the rawdir string
    size_t length = strlen(rawdir);
    // Allocate memory for the newdir string
    char *newdir = malloc((length + 1) * sizeof(char));
    if (newdir == NULL) {
        perror("Memory allocation failed");
        return NULL;
    }
    // Copy the rawdir string into newdir
    strcpy(newdir, rawdir);
    // Perform the conversion based on the operating system
    for (size_t i = 0; i < length; i++) {
        #ifdef _WIN32
            if (newdir[i] == '/') {
                newdir[i] = '\\';
            }
        #else
            if (newdir[i] == '\\') {
                newdir[i] = '/';
            }
        #endif
    }
    // The user is responsible to free (newdir) after function call
    return newdir;
}

void OSmkdir(char *path) {
    #ifdef _WIN32
        if (_mkdir(path) != 0) {
            // If error is different than the existance of the directory, return..
            if (errno != EEXIST) {
                perror(path);   // Print error
                return;
            }
        }  
    #else
        if (mkdir(path, S_IRWXU) != 0) {
            // If error is different than the existance of the directory, return..
            if (errno != EEXIST) {
                perror(path);   // Print error
                return;
            }
        }
    #endif
}

void create_dir(const char *pathname) {
    errno = 0;
    // If pathname doesn't exists, return error
    if (pathname == NULL) {
        perror(pathname);
        return;
    }
    // Declare variables
    const size_t len = strlen(pathname);                    // Size of the pathname
    char *path = malloc(MAX_FILENAME_LEN * sizeof(char));   // Allocate temporary memory to copy pathname
    char *ch;                                               // Pointer to each character in the string
    // Check if len is too long to be stored into path
    if (len > MAX_FILENAME_LEN * sizeof(char) - 1) {
        errno = ENAMETOOLONG;                               // Return error if its too long
        perror(pathname);                                   // Print error
        return;
    }
    // Copy string pathname into string path
    strcpy(path, pathname);
    // Read the string char by char
    for (ch = path + 1; *ch; ch++) {
        // Check if ch is a separator
        if (*ch == SEP) {
            // Truncate temporarily the string to check if dir exists
            *ch = '\0';
            // Check if exists. If not, it is created
            OSmkdir(path);
            // Insert the separator again in the string to continue
            *ch = SEP;
        }
    }
    
    // Check again the full path to be sure. If dont exist, create
    OSmkdir(path);   
    free(path);
}

void directory_exists(const char *pathname) {
    // Check if the directory exists
    if (stat(pathname, &info) != 0) {
        printf("   Cannot find output directory\n   Creating output directory...\n");
        create_dir(pathname);
    } 
}

bool file_exists(const char* filename) {
    // Try to open file with same name as filename
    FILE *testfile = fopen(filename, "r");
    // Check the existance of a file called filename
    if (testfile != NULL) {
        // Close testfile and returns true if the file exists
        fclose(testfile); 
        return true;
    }
    else {
        // Reurn false if the file does not exist
        return false;
    }
}

FILE *create_output_file(char *name, const char *ext, const char *dir) {
    // Check for safety
    if ((name == NULL) || (ext == NULL) || (dir == NULL)) {
        print_error("Failed to create output file.\n");
        print_error("Please check the code.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Check if directory (dir) exists. If not, create
    directory_exists(dir);
    // Creates a variable storing the total size of name (200) 
    size_t size = sizeof(name);
    // Creates a variable storing the size of the first attributed name
    size_t len = strlen(name);
    // Insert the extension of the file in the end of the name string
    snprintf(name + len, size - len, "%s", ext);
    // Checks if a file with the same name already exists in the directory
    int j = 0;
    while (file_exists(name)) {
        j++;
        if (j >= MAX_FILE_NUM) {
            // All possible file versions in the limit of MAX_FILE_NUM tried
            return NULL;
        }
        snprintf(name + len, size - len, "(%i)%s", j, ext);
    } 

    return fopen(name, "w");
}

FILE *name_and_create_output_files(const char *systemname, const char *directory, const char *module, const char *ext) {
    // Create output files to store results
    char output_filename[MAX_FILENAME_LEN + 1];
    // Convert directory string to match operational system
    char *dir = convert_dir(directory);
    snprintf(output_filename, sizeof(output_filename), "%s%s_%s", dir, systemname, module); // Assign name for output file without extension
    FILE *output_file = create_output_file(output_filename, ext, dir);                       // Create output file 
    // Free memory
    free(dir);

    return output_file;
}

char* get_input_filename(void) {
    char* filename = malloc((MAX_FILENAME_LEN + 1) * sizeof(char));
    if (filename == NULL) {
        print_error("Failed to allocate memory for the input filename.\n");
        return NULL;
    }
    // Clear the input buffer
    int c;
    while ((c = getchar()) != '\n' && c != EOF);

    // Asks for the user for the filename
    printf("Enter Input Filename: ");
    fflush(stdout);  // Flush stdout to ensure prompt is displayed
    // Read input using fgets
    if (fgets(filename, MAX_FILENAME_LEN + 1, stdin) == NULL) {
        print_error("Failed to read input.\n");
        free(filename);
        return NULL;
    }
    // Remove trailing newline character, if present
    size_t length = strlen(filename);
    if (length > 0 && filename[length - 1] == '\n') {
        filename[length - 1] = '\0';
        length--;
    }
    // Check if input exceeds maximum length
    if (length == MAX_FILENAME_LEN && getchar() != '\n') {
        print_error("Invalid input: Filename exceeds maximum length of %d characters.\n", MAX_FILENAME_LEN);
        print_error("Please, reduce the name of the file or/and allocate the file in a directory with smaller string size\n");
        free(filename);
        return NULL;
    }

    return filename;
}


// Write Results for Nonlinear Dynamics Toolbox

void write_results(FILE *output_file, int dim, double t, double *x, int mode) {
    // Check the mode of the function
    if (mode == 0) {
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 1) {
        // Header
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "\n");
        // Initial Conditions
        fprintf(output_file, "%.10f ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        fprintf(output_file, "\n");
    } 
    else if (mode == 2) {
        fprintf(output_file, "%.10f ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        fprintf(output_file, "\n");
    }
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void write_results_lyap(FILE *output_file, int dim, double t, double *lambda, double *s_lambda, int mode) {
    // Check the mode of the function
    if (mode == 1) {
        // Header
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%i] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "sLE[%i] ", i);
        }
        fprintf(output_file, "\n");
        // Initial Exponents
        fprintf(output_file, "%.10f ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", s_lambda[i]);
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 2) {
        fprintf(output_file, "%.10f ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", s_lambda[i]);
        }
        fprintf(output_file, "\n");
    } 
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void write_results_ftimeseries(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int mode) {
    // Check the mode of the function
    if (mode == 1) {
        // Header
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%i] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "sLE[%i] ", i);
        }
        fprintf(output_file, "\n");
        // Initial Conditions
        fprintf(output_file, "%.10lf ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", s_lambda[i]);
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 2) {
        fprintf(output_file, "%.10lf ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", s_lambda[i]);
        }
        fprintf(output_file, "\n");
    } 
    else if (mode == 3) {
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 4) {
        fprintf(output_file, "%.10lf ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        fprintf(output_file, "\n");
    }
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void write_bifurc_results(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, int mode) {
    // Check the mode of the function
    if (mode == 0) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "\n");
    } 
    else if (mode == 1) {
        fprintf(output_file, "%.10f ", varpar);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 2) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "xmax[%i] ", i);
            fprintf(output_file, "xmin[%i] ", i);
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 3) {
        fprintf(output_file, "%.10f ", varpar);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", xmax[i]);
            fprintf(output_file, "%.10lf ", xmin[i]);
        }
        fprintf(output_file, "\n");
    }
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void write_fbifurc_results(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *LE, int attractor, double **poinc, int diffattrac, int mode) {
    // Check the mode of the function
    // Header for poincare bifurc
    if (mode == 0) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "Attractor DiffAttrac\n");
    } 
    // Write Results for poincare bifurc
    else if (mode == 1) {
        for(int q = 0; q < np - trans; q ++) {
            fprintf(output_file, "%.10lf ", varpar);
            for (int i = 0; i < dim; i++) {
                fprintf(output_file, "%.10lf ", poinc[q][i]);
            }
            fprintf(output_file, "%d %d\n", attractor, diffattrac);
        }
    }
    // Header for lyapunov, min, max values
    else if (mode == 2) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "xmin[%d] ", i);
            fprintf(output_file, "xmax[%d] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%d] ", i);
        }
        fprintf(output_file, "Attractor DiffAttrac\n");
    }
    else if (mode == 3) {
        fprintf(output_file, "%.10f ", varpar);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", xmax[i]);
            fprintf(output_file, "%.10lf ", xmin[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", LE[i]);
        }
        fprintf(output_file, "%d %d\n", attractor, diffattrac);
    }
    // Error
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void write_dyndiag_results(FILE *output_file, int dim, double varparX, double varparY, int attractor, double *LE, int diffattrac, int mode) {
    // Check the mode of the function
    // Header
    if (mode == 0) {
        fprintf(output_file, "%s %s %s %s ", "CparY", "CparX", "Attractor", "DiffAttrac");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%i] ", i);
        }
        fprintf(output_file, "\n");
    } 
    // Write Results
    else if (mode == 1) {
        fprintf(output_file, "%.10lf %.10lf %i %i ", varparY, varparX, attractor, diffattrac);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", LE[i]);
        }
        fprintf(output_file, "\n");
    }
    // Error
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void p_write_fdyndiag_results(FILE *output_file, int dim, double **results, int pixels) {
    // Header
    fprintf(output_file, "%s %s %s %s ", "CparY", "CparX", "Attractor", "DiffAttrac");
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "LE[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "xmax[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "xmin[%d] ", i);
    }
    fprintf(output_file, "\n");
    // Write Results
    for (int i = 0; i < pixels; i++) {
        fprintf(output_file, "%.10lf %.10lf %d %d ", results[i][0], results[i][1], (int)results[i][2], (int)results[i][3]);
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4]);
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+dim]);
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+(2*dim)]);
        }
        fprintf(output_file, "\n");
    }
}

void p_write_dyndiag_results(FILE *output_file, int dim, double **results, int pixels) {
    // Header
    fprintf(output_file, "%s %s %s %s ", "CparY", "CparX", "Attractor", "DiffAttrac");
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "xmax[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "xmin[%d] ", i);
    }
    fprintf(output_file, "\n");
    // Write Results
    for (int i = 0; i < pixels; i++) {
        fprintf(output_file, "%.10lf %.10lf %d %d ", results[i][0], results[i][1], (int)results[i][2], (int)results[i][3]);
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4]);
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+dim]);
        }
        fprintf(output_file, "\n");
    }
}

void p_write_epbasin_results(FILE *output_file, double **results, int pixels, int dim) {
    printf("\n  Writing Results in Output File...\n");
    // Header
    fprintf(output_file, "%s %s ", "CparY", "CparX");
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "x[%d] ", i);
    }
    fprintf(output_file, "%s\n", "Attractor");
    // Write Results
    for (int i = 0; i < pixels; i++) {
        fprintf(output_file, "%.10lf %.10lf ", results[i][0], results[i][1]);
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+2]);
        }
        fprintf(output_file, "%d\n", (int)results[i][dim + 2]);
    }
}

// Write results for Energy Harvesting Toolbox

void EH_write_bifurc_results(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int nrms, int *rmsindex, double *xrms, double *overallxrms,
                             int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, int mode) {
    // Check the mode of the function
    if (mode == 0) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "\n");
    } 
    else if (mode == 1) {
        fprintf(output_file, "%.10f ", varpar);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 2) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "xMIN[%d] ", i);
            fprintf(output_file, "xMAX[%d] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "OverallxMIN[%d] ", i);
            fprintf(output_file, "OverallxMAX[%d] ", i);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "xRMS[%d] ", rmsindex[i]);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "OverallxRMS[%d] ", rmsindex[i]);
        }
        // If it has any custom calculations, add custom header
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%s ", customnames[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 3) {
        fprintf(output_file, "%.10f ", varpar);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", xmin[i]);
            fprintf(output_file, "%.10lf ", xmax[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", overallxmin[i]);
            fprintf(output_file, "%.10lf ", overallxmax[i]);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "%.10lf ", xrms[rmsindex[i]]);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "%.10lf ", overallxrms[rmsindex[i]]);
        }
        // If it has any custom calculations, add custom values
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%.10lf ", customvalue[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
    }
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void EH_write_fbifurc_results(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *overallxmin, double *overallxmax, double *LE, int attractor, double **poinc, int diffattrac, int nrms, int *rmsindex, double *xrms, double *overallxrms,
                              int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, int mode) {
    // Check the mode of the function
    // Header for poincare bifurc
    if (mode == 0) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "Attractor DiffAttrac\n");
    } 
    // Write Results for poincare bifurc
    else if (mode == 1) {
        for(int q = 0; q < np - trans; q ++) {
            fprintf(output_file, "%.10lf ", varpar);
            for (int i = 0; i < dim; i++) {
                fprintf(output_file, "%.10lf ", poinc[q][i]);
            }
            fprintf(output_file, "%d %d\n", attractor, diffattrac);
        }
    }
    // Header for lyapunov, min, max values
    else if (mode == 2) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "xMIN[%d] ", i);
            fprintf(output_file, "xMAX[%d] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "OverallxMIN[%d] ", i);
            fprintf(output_file, "OverallxMAX[%d] ", i);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "xRMS[%d] ", rmsindex[i]);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "OverallxRMS[%d] ", rmsindex[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%d] ", i);
        }
        fprintf(output_file, "Attractor DiffAttrac ");
        // If it has any custom calculations, add custom header
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%s ", customnames[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 3) {
        fprintf(output_file, "%.10f ", varpar);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", xmin[i]);
            fprintf(output_file, "%.10lf ", xmax[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", overallxmin[i]);
            fprintf(output_file, "%.10lf ", overallxmax[i]);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "%.10lf ", xrms[rmsindex[i]]);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "%.10lf ", overallxrms[rmsindex[i]]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", LE[i]);
        }
        fprintf(output_file, "%d %d ", attractor, diffattrac);
        // If it has any custom calculations, add custom values
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%.10lf ", customvalue[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
    }
    // Error
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void EH_p_write_fdyndiag_results(FILE *output_file, int dim, int nrms, int *rmsindex, double **results, int pixels, int ncustomvalues, char **customnames, int nprintf, int *printfindex) {
    // Header
    fprintf(output_file, "%s %s %s %s ", "CparY", "CparX", "Attractor", "DiffAttrac");
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "LE[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "xMAX[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "xMIN[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "OverallxMAX[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "OverallxMIN[%d] ", i);
    }
    if (nrms > 0) {
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "xRMS[%d] ", rmsindex[i]);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "OverallxRMS[%d] ", rmsindex[i]);
        }
    }
    // If it has any custom calculations, add custom header
    if (ncustomvalues > 0) { 
        for (int i = 0; i < nprintf; i++) {
            fprintf(output_file, "%s ", customnames[printfindex[i]]);
        } 
    }
    fprintf(output_file, "\n");
    // Write Results
    for (int i = 0; i < pixels; i++) {
        fprintf(output_file, "%.10lf %.10lf %d %d ", results[i][0], results[i][1], (int)results[i][2], (int)results[i][3]);
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4]); // LE
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+dim]); // xMAX
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+(2*dim)]); // xMin
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+(3*dim)]); // overallxmax
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+(4*dim)]); // overallxmin
        }
        if (nrms > 0) {
            for (int j = 0; j < nrms; j++) {
                fprintf(output_file, "%.10lf ", results[i][j+4+(5*dim)]); // xRMS
            }
            for (int j = 0; j < nrms; j++) {
                fprintf(output_file, "%.10lf ", results[i][j+4+(5*dim)+nrms]); //OverallxRMS
            }
        }
        // If it has any custom calculations, add custom values
        if (ncustomvalues > 0) {
            for (int j = 0; j < nprintf; j++) {
                fprintf(output_file, "%.10lf ", results[i][j+4+(5*dim)+(2*nrms)]);  //customvalues
            }
        }
        fprintf(output_file, "\n");
    }
}

void EH_p_write_dyndiag_results(FILE *output_file, int dim, int nrms, int *rmsindex, double **results, int pixels, int ncustomvalues, char **customnames, int nprintf, int *printfindex) {
    // Header
    fprintf(output_file, "%s %s %s %s ", "CparY", "CparX", "Attractor", "DiffAttrac");
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "xMAX[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "xMIN[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "OverallxMAX[%d] ", i);
    }
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "OverallxMIN[%d] ", i);
    }
    if (nrms > 0) {
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "xRMS[%d] ", rmsindex[i]);
        }
        for (int i = 0; i < nrms; i++) {
            fprintf(output_file, "OverallxRMS[%d] ", rmsindex[i]);
        }
    }
    // If it has any custom calculations, add custom header
    if (ncustomvalues > 0) { 
        for (int i = 0; i < nprintf; i++) {
            fprintf(output_file, "%s ", customnames[printfindex[i]]);
        } 
    }
    fprintf(output_file, "\n");
    // Write Results
    for (int i = 0; i < pixels; i++) {
        fprintf(output_file, "%.10lf %.10lf %d %d ", results[i][0], results[i][1], (int)results[i][2], (int)results[i][3]);
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4]); // xmax
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+dim]); // xmin
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+(2*dim)]); // overallxmax
        }
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4+(3*dim)]); // overallxmin
        }
        if (nrms > 0) {
            for (int j = 0; j < nrms; j++) {
                fprintf(output_file, "%.10lf ", results[i][j+4+(4*dim)]);
            }
            for (int j = 0; j < nrms; j++) {
                fprintf(output_file, "%.10lf ", results[i][j+4+(4*dim)+nrms]);
            }
        }
        // If it has any custom calculations, add custom values
        if (ncustomvalues > 0) {
            for (int j = 0; j < nprintf; j++) {
                fprintf(output_file, "%.10lf ", results[i][j+4+(4*dim)+(2*nrms)]);
            }
        }
        fprintf(output_file, "\n");
    }
}

void EH_write_timeseries_results(FILE *output_file, int dim, double t, double *x, int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, int mode) {
    // Check the mode of the function
    if (mode == 1) {
        // Add Header
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        // If it has any custom calculations, add custom header
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%s ", customnames[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
        // Add Initial Conditions
        fprintf(output_file, "%.10f ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        // If it has any custom calculations, add custom values
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%.10lf ", customvalue[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
    } 
    else if (mode == 2) {
        // Add results
        fprintf(output_file, "%.10f ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        // If it has any custom calculations, add custom results
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%.10lf ", customvalue[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
    }
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void EH_write_ftimeseries_results(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, int mode) {
    // Check the mode of the function
    if (mode == 1) {
        // Header
        fprintf(output_file, "Time ");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%i] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "sLE[%i] ", i);
        }
        // If it has any custom calculations, add custom header
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%s ", customnames[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
        // Initial Conditions
        fprintf(output_file, "%.10lf ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", s_lambda[i]);
        }
        // If it has any custom calculations, add custom values
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%.10lf ", customvalue[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
    }
    else if (mode == 2) {
        // Add Results
        fprintf(output_file, "%.10lf ", t);
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", s_lambda[i]);
        }
        // If it has any custom calculations, add custom values
        if (ncustomvalues > 0) { 
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%.10lf ", customvalue[printfindex[i]]);
            } 
        }
        fprintf(output_file, "\n");
    } 
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}