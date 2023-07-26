#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <math.h>
#include <sys/stat.h>
#include "defines.h"
#include "msg.h"

#ifdef _WIN32
    #include <direct.h>
    const char SEP = '\\';
#else
    const char SEP = '/';
#endif    

struct stat info;

// Functions to handle directories, file creation, etc

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

// Functions to write specific blocks of information into output files
static void write_time(FILE *output_file, double t, int mode) {
    // Write header
    if (mode == 1) {
        fprintf(output_file, "Time ");
    }
    // Write results
    else if (mode == 2) {
        fprintf(output_file, "%.10f ", t);
    }
    else {
        print_debug("Failed to write results in output file with function 'write_time()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_state_vars(FILE *output_file, int dim, double *x, int mode) {
    // Write header
    if (mode == 1) {
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%d] ", i);
        }
    }
    // Write results
    else if (mode == 2) {
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", x[i]);
        }
    }
    else {
        print_debug("Failed to write results in output file with function 'write_state_vars()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_lyap_exp(FILE *output_file, int dim, double *lambda, double *s_lambda, bool short_lamb, int mode) {
    // Write header
    if (mode == 1) {
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%d] ", i);
        }
        if (short_lamb == true) {
            for (int i = 0; i < dim; i++) {
                fprintf(output_file, "sLE[%d] ", i);
            }
        }
    }
    // Write results
    else if (mode == 2) {
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", lambda[i]);
        }
        if (short_lamb == true) {
            for (int i = 0; i < dim; i++) {
                fprintf(output_file, "%.10lf ", s_lambda[i]);
            }
        }
    }
    else {
        print_debug("Failed to write results in output file with function 'write_lyap_exp()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_state_vars_angles(FILE *output_file, double *x, ang_info *angles, int mode) {
    if (angles->n_angles > 0) {
        // Write header
        if (mode == 1) {
            for(int i = 0; i < angles->n_angles; i++) {
                fprintf(output_file, "x[%d]_remainder ", angles->index[i]);
            }
        }
        // Write results
        else if (mode == 2) {
            for(int i = 0; i < angles->n_angles; i++) {
                fprintf(output_file, "%.10lf ", remainder(x[angles->index[i]], TWOPI));
            }
        }
        else {
            print_debug("Failed to write results in output file with function 'write_state_vars()' using mode (%d)...\n", mode);
            return;
        }
    }
    else {
        return;
    }
}

static void write_min_max(FILE *output_file, int dim, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int mode) {
    // Header
    if (mode == 1) {
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "xMIN[%d] ", i);
            fprintf(output_file, "xMAX[%d] ", i);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "OverallxMIN[%d] ", i);
            fprintf(output_file, "OverallxMAX[%d] ", i);
        }
    } 
    // Results
    else if (mode == 2) {
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", xmin[i]);
            fprintf(output_file, "%.10lf ", xmax[i]);
        }
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "%.10lf ", overallxmin[i]);
            fprintf(output_file, "%.10lf ", overallxmax[i]);
        }
    }
    else {
        print_debug("Failed to write results in output file with function 'write_min_max()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_min_max_angles(FILE *output_file, ang_info *angles, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int mode) {
    // If it has any angle within the set of state variables, add remainder (remainder of value/2*pi)
    if (angles->n_angles > 0) {
         // Header results
        if (mode == 1) {
            for(int i = 0; i < angles->n_angles; i++) {
                fprintf(output_file, "xMIN[%d]_remainder ", angles->index[i]);
                fprintf(output_file, "xMAX[%d]_remainder ", angles->index[i]);
            }
            for(int i = 0; i < angles->n_angles; i++) {
                fprintf(output_file, "OverallxMIN[%d]_remainder ", angles->index[i]);
                fprintf(output_file, "OverallxMAX[%d]_remainder ", angles->index[i]);
            }
        }
        // Results 
        else if (mode == 2) {
            // xmin and xmax
            for (int i = 0; i < angles->n_angles; i++) {
                if (fabs(xmax[angles->index[i]] - xmin[angles->index[i]]) > TWOPI) {
                    fprintf(output_file, "%.10lf ", -PI);
                    fprintf(output_file, "%.10lf ", PI);
                }
                else {
                    fprintf(output_file, "%.10lf ", remainder(xmin[angles->index[i]], TWOPI));
                    fprintf(output_file, "%.10lf ", remainder(xmax[angles->index[i]], TWOPI));
                }
            }
            // Overallxmin and Overallxmax
            for (int i = 0; i < angles->n_angles; i++) {
                if (fabs(overallxmax[angles->index[i]] - overallxmin[angles->index[i]]) > TWOPI) {
                    fprintf(output_file, "%.10lf ", -PI);
                    fprintf(output_file, "%.10lf ", PI);
                }
                else {
                    fprintf(output_file, "%.10lf ", remainder(overallxmin[angles->index[i]], TWOPI));
                    fprintf(output_file, "%.10lf ", remainder(overallxmax[angles->index[i]], TWOPI));
                }
            }
        }
        else {
            print_debug("Failed to write results in output file with function 'write_max_min_angles()' using mode (%d)...\n", mode);
            return;
        }
    }
    else {
        return;
    }
}

static void write_rms_state_vars(FILE *output_file, int nrms, double *xrms, double *overallxrms, int *rmsindex, int mode) {
    // If there is any RMS calculations to be performed
    if (nrms > 0) {
        // Header
        if (mode == 1) {
            for (int i = 0; i < nrms; i++) {
                fprintf(output_file, "xRMS[%d] ", rmsindex[i]);
            }
            for (int i = 0; i < nrms; i++) {
                fprintf(output_file, "OverallxRMS[%d] ", rmsindex[i]);
            }
        }
        // Results
        else if (mode == 2) {
            for (int i = 0; i < nrms; i++) {
                fprintf(output_file, "%.10lf ", xrms[rmsindex[i]]);
            }
            for (int i = 0; i < nrms; i++) {
                fprintf(output_file, "%.10lf ", overallxrms[rmsindex[i]]);
            }
        }
        else {
            print_debug("Failed to write results in output file with function 'write_rms_state_vars()' using mode (%d)...\n", mode);
            return;
        }
    }
    else {
        return;
    }
}

static void write_customcalc(FILE *output_file, int ncustomvalues, double *customvalue, int nprintf, char **customnames, int *printfindex, int mode) {
    // If it has any custom calculations, add custom values
    if (ncustomvalues > 0) { 
        // Header
        if (mode == 1) {
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%s ", customnames[printfindex[i]]);
            } 
        }
        // Results
        else if (mode == 2) {
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%.10lf ", customvalue[printfindex[i]]);
            } 
        }
        else {
            print_debug("Failed to write results in output file with function 'write_customcalc()' using mode (%d)...\n", mode);
            return;
        }
    }
    else {
        return;
    }
}

static void write_bifurc_control_parameter(FILE *output_file, double varpar, int mode) {
    // Write header
    if (mode == 1) {
        fprintf(output_file, "Cpar ");
    }
    // Write results
    else if (mode == 2) {
        fprintf(output_file, "%.10f ", varpar);
    }
    else {
        print_debug("Failed to write results in output file with function 'write_bifurc_control_parameter()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_attractor(FILE *output_file, int attractor, int mode) {
    // Header
    if (mode == 1) {
        fprintf(output_file, "%s ", "Attractor");
    } 
    else if (mode == 2) {
        fprintf(output_file, "%d ", attractor);
    }
    else {
        print_debug("Failed to write results in output file with function 'write_attractor()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_2D_control_params_result_matrix_in_file(FILE *output_file, int *col_offset, double *results, int mode) {
    // Header
    if (mode == 1) {
        fprintf(output_file, "CparY CparX ");
    } 
    else if (mode == 2) {
        fprintf(output_file, "%.10lf %.10lf ", results[0], results[1]);
        (*col_offset) += 2;
    }
    else {
        print_debug("Failed to write results in output file with function 'write_2D_control_params_result_matrix_in_file()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_attractor_results_matrix_in_file(FILE *output_file, int *col_offset, double *results, int mode) {
    // Header
    if (mode == 1) {
        fprintf(output_file, "%s ", "Attractor");
    } 
    else if (mode == 2) {
            fprintf(output_file, "%d ", (int)results[2]);
            (*col_offset) += 1;
    }
    else {
        print_debug("Failed to write results in output file with function 'write_attractor_results_matrix_in_file()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_min_max_result_matrix_in_file(FILE *output_file, int dim, int *col_offset, double *results, int mode) {
    // Header
    if (mode == 1) {
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
    }
    // Results
    else if (mode == 2) {
        // Write xmax[dim]
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[(*col_offset) + j]);
        }
        (*col_offset) += dim;
        // Write xmin[dim]
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[(*col_offset) + j]); 
        }
        (*col_offset) += dim;
        // Write overallxmax[dim]
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[(*col_offset) + j]); 
        }
        (*col_offset) += dim;
        // Write overallxmin[dim]
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[(*col_offset) + j]); 
        }
        (*col_offset) += dim;
    }
    else {
        print_debug("Failed to write results in output file with function 'write_min_max_result_matrix_in_file()' using mode (%d)...\n", mode);
        return;
    }
}

static void write_rms_result_matrix_in_file(FILE *output_file, int nrms, int *rmsindex, int *col_offset, double *results, int mode) {
    if (nrms > 0) {
        // Header
        if (mode == 1) {
            for (int i = 0; i < nrms; i++) {
                fprintf(output_file, "xRMS[%d] ", rmsindex[i]);
            }
            for (int i = 0; i < nrms; i++) {
                fprintf(output_file, "OverallxRMS[%d] ", rmsindex[i]);
            }
        }
        // Results
        else if (mode == 2) {
            // Write xRMS
            for (int j = 0; j < nrms; j++) {
                fprintf(output_file, "%.10lf ", results[(*col_offset) + j]);
            }
            (*col_offset) += nrms;
            // Write OverallxRMS
            for (int j = 0; j < nrms; j++) {
                fprintf(output_file, "%.10lf ", results[(*col_offset) + j]);
            }
            (*col_offset) += nrms;
        }
        else {
            print_debug("Failed to write results in output file with function 'write_rms_result_matrix_in_file()' using mode (%d)...\n", mode);
            return;
        }
    }
    else {
        return;
    }
}

static void write_customcalc_result_matrix_in_file(FILE *output_file, int ncustomvalues, int nprintf, int *printfindex, char **customnames, int *col_offset, double *results, int mode) {
    if (ncustomvalues > 0) {
        // Header
        if (mode == 1) {
            for (int i = 0; i < nprintf; i++) {
                fprintf(output_file, "%s ", customnames[printfindex[i]]);
            }
        }
        // Results
        else if (mode == 2) {
            for (int j = 0; j < nprintf; j++) {
                fprintf(output_file, "%.10lf ", results[(*col_offset) + j]);
            }
            (*col_offset) += nprintf;
        }
        else {
            print_debug("Failed to write results in output file with function 'write_customcalc_result_matrix_in_file()' using mode (%d)...\n", mode);
            return;
        }
    }
    else {
        return;
    }
}

static void write_min_max_angles_result_matrix_in_file(FILE *output_file, ang_info *angles, int *col_offset, double *results, int mode) {
    // Header
    if (angles->n_angles > 0)  {
        if (mode == 1) {
                for (int i = 0; i < angles->n_angles; i++) {
                    fprintf(output_file, "xMAX[%d]_remainder ", angles->index[i]);
                }
                for (int i = 0; i < angles->n_angles; i++) {
                    fprintf(output_file, "xMIN[%d]_remainder ", angles->index[i]);
                }
                for (int i = 0; i < angles->n_angles; i++) {
                    fprintf(output_file, "OverallxMAX[%d]_remainder ", angles->index[i]);
                }
                for (int i = 0; i < angles->n_angles; i++) {
                    fprintf(output_file, "OverallxMIN[%d]_remainder ", angles->index[i]);
                }
            }
            // Results
            else if (mode == 2) {
                // Write xmax[angles->n_angles]_remainder
                for (int j = 0; j < angles->n_angles; j++) {
                    fprintf(output_file, "%.10lf ", results[(*col_offset) + j]);
                }
                (*col_offset) += angles->n_angles;
                // Write xmin[angles->n_angles]_remainder
                for (int j = 0; j < angles->n_angles; j++) {
                    fprintf(output_file, "%.10lf ", results[(*col_offset) + j]); 
                }
                (*col_offset) += angles->n_angles;
                // Write overallxmax[angles->n_angles]_remainder
                for (int j = 0; j < angles->n_angles; j++) {
                    fprintf(output_file, "%.10lf ", results[(*col_offset) + j]); 
                }
                (*col_offset) += angles->n_angles;
                // Write overallxmin[angles->n_angles]_remainder
                for (int j = 0; j < angles->n_angles; j++) {
                    fprintf(output_file, "%.10lf ", results[(*col_offset) + j]); 
                }
                (*col_offset) += angles->n_angles;
            }
            else {
                print_debug("Failed to write results in output file with function 'write_min_max_angles_result_matrix_in_file()' using mode (%d)...\n", mode);
                return;
            }
    }
    else {
        return;
    }
}

static void write_LE_result_matrix_in_file(FILE *output_file, int dim, int *col_offset, double *results, int mode) {
    // Header
    if (mode == 1) {
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "LE[%d] ", i);
        }
    }
    // Results
    else if (mode == 2) {
        // Write LE[dim]
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[(*col_offset) + j]);
        }
        (*col_offset) += dim;
    }
    else {
        print_debug("Failed to write results in output file with function 'write_LE_result_matrix_in_file()' using mode (%d)...\n", mode);
        return;
    }
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

// Write results for Harmonic Oscillator (HOS) Toolbox

void HOS_write_timeseries_results(FILE *output_file, int dim, double t, double *x, int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, ang_info *angles, int mode) {
    // Check the mode of the function
    if (mode == 1) {
        // Add Header
        write_time(output_file, t, 1);
        write_state_vars(output_file, dim, x, 1);
        write_state_vars_angles(output_file, x, angles, 1);
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 1);
        fprintf(output_file, "\n");
        // Add Initial Conditions
        write_time(output_file, t, 2);
        write_state_vars(output_file, dim, x, 2);
        write_state_vars_angles(output_file, x, angles, 2);
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 2);
        fprintf(output_file, "\n");
    } 
    else if (mode == 2) {
        write_time(output_file, t, 2);
        write_state_vars(output_file, dim, x, 2);
        write_state_vars_angles(output_file, x, angles, 2);
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 2);
        fprintf(output_file, "\n");
    }
    else {
        print_debug("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void HOS_write_poinc_results(FILE *output_file, int dim, double t, double *x, ang_info *angles, int mode) {
    // Check the mode of the function
    if (mode == 1) {
        // Add Header
        write_time(output_file, t, 1);
        write_state_vars(output_file, dim, x, 1);
        write_state_vars_angles(output_file, x, angles, 1);
        fprintf(output_file, "\n");
    } 
    else if (mode == 2) {
        // Add results
        write_time(output_file, t, 2);
        write_state_vars(output_file, dim, x, 2);
        write_state_vars_angles(output_file, x, angles, 2);
        fprintf(output_file, "\n");
    }
    else {
        print_debug("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void HOS_write_ftimeseries_results(FILE *output_file, int dim, double t, double *x, double *lambda, double *s_lambda, int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, ang_info *angles, int mode) {
    // Check the mode of the function
    if (mode == 1) {
        // Header
        write_time(output_file, t, 1);
        write_state_vars(output_file, dim, x, 1);
        write_state_vars_angles(output_file, x, angles, 1);
        write_lyap_exp(output_file, dim, lambda, s_lambda, true, 1);
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 1);
        fprintf(output_file, "\n");
        // Initial Conditions
        write_time(output_file, t, 2);
        write_state_vars(output_file, dim, x, 2);
        write_state_vars_angles(output_file, x, angles, 2);
        write_lyap_exp(output_file, dim, lambda, s_lambda, true, 2);
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 2);
        fprintf(output_file, "\n");
    }
    else if (mode == 2) {
        write_time(output_file, t, 2);
        write_state_vars(output_file, dim, x, 2);
        write_state_vars_angles(output_file, x, angles, 2);
        write_lyap_exp(output_file, dim, lambda, s_lambda, true, 2);
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 2);
        fprintf(output_file, "\n");
    } 
    else {
        print_debug("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void HOS_write_bifurc_results(FILE *output_file, int dim, double varpar, double *x, double *xmin, double *xmax, double *overallxmin, double *overallxmax, int nrms, int *rmsindex, double *xrms, double *overallxrms,
                             int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, ang_info *angles, int mode) {
    // Header of the poincare file
    if (mode == 0) {
        write_bifurc_control_parameter(output_file, varpar, 1);
        write_state_vars(output_file, dim, x, 1);
        write_state_vars_angles(output_file, x, angles, 1);
        fprintf(output_file, "\n");
    } 
    // Results of the poincare file
    else if (mode == 1) {
        write_bifurc_control_parameter(output_file, varpar, 2);
        write_state_vars(output_file, dim, x, 2);
        write_state_vars_angles(output_file, x, angles, 2);
        fprintf(output_file, "\n");
    }
    // Header of measurements file
    else if (mode == 2) {
        write_bifurc_control_parameter(output_file, varpar, 1);
        write_min_max(output_file, dim, xmin, xmax, overallxmin, overallxmax, 1);
        write_min_max_angles(output_file, angles, xmin, xmax, overallxmin, overallxmax, 1);
        write_rms_state_vars(output_file, nrms, xrms, overallxrms, rmsindex, 1);
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 1);
        fprintf(output_file, "\n");
    }
    // Results of measurements file
    else if (mode == 3) {
        write_bifurc_control_parameter(output_file, varpar, 2);
        write_min_max(output_file, dim, xmin, xmax, overallxmin, overallxmax, 2);
        write_min_max_angles(output_file, angles, xmin, xmax, overallxmin, overallxmax, 2);
        write_rms_state_vars(output_file, nrms, xrms, overallxrms, rmsindex, 2);
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 2);
        fprintf(output_file, "\n");
    }
    else {
        print_debug("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void HOS_write_fbifurc_results(FILE *output_file, int dim, int np, int trans, double varpar, double *x, double *xmin, double *xmax, double *overallxmin, double *overallxmax, double *LE, int attractor, size_t npoinc, double **poinc, int nrms, int *rmsindex, double *xrms, double *overallxrms,
                              int ncustomvalues, char **customnames, double *customvalue, int nprintf, int *printfindex, ang_info *angles, int mode) {
    // Check the mode of the function
    // Header for poincare bifurc
    if (mode == 0) {
        //write_fbifurc_poinc(output_file, dim, varpar, npoinc, poinc, attractor, angles, 1);
        write_bifurc_control_parameter(output_file, varpar, 1);
        write_state_vars(output_file, dim, NULL, 1);
        write_state_vars_angles(output_file, NULL, angles, 1);
        write_attractor(output_file, attractor, 1);
        fprintf(output_file, "\n");
    } 
    // Write Results for poincare bifurc
    else if (mode == 1) {
        for(int q = 0; q < npoinc; q ++) {
            //write_fbifurc_poinc(output_file, dim, varpar, npoinc, poinc, attractor, angles, 2);
            write_bifurc_control_parameter(output_file, varpar, 2);
            write_state_vars(output_file, dim, poinc[q], 2);
            write_state_vars_angles(output_file, poinc[q], angles, 2);
            write_attractor(output_file, attractor, 2);
            fprintf(output_file, "\n");
        }
    }
    // Header for lyapunov, min, max values, customvalues
    else if (mode == 2) {
        write_bifurc_control_parameter(output_file, varpar, 1);
        write_min_max(output_file, dim, xmin, xmax, overallxmin, overallxmax, 1);
        write_min_max_angles(output_file, angles, xmin, xmax, overallxmin, overallxmax, 1);
        write_rms_state_vars(output_file, nrms, xrms, overallxrms, rmsindex, 1);
        write_lyap_exp(output_file, dim, LE, NULL, false, 1);
        write_attractor(output_file, attractor, 1);        
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 1);
        fprintf(output_file, "\n");
    }
    // Write Results for lyapunov, min, max values, customvalues
    else if (mode == 3) {
        write_bifurc_control_parameter(output_file, varpar, 2);
        write_min_max(output_file, dim, xmin, xmax, overallxmin, overallxmax, 2);
        write_min_max_angles(output_file, angles, xmin, xmax, overallxmin, overallxmax, 2);
        write_rms_state_vars(output_file, nrms, xrms, overallxrms, rmsindex, 2);
        write_lyap_exp(output_file, dim, LE, NULL, false, 2);
        write_attractor(output_file, attractor, 2);        
        write_customcalc(output_file, ncustomvalues, customvalue, nprintf, customnames, printfindex, 2);
        fprintf(output_file, "\n");
    }
    // Error
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void HOS_p_write_dyndiag_results(FILE *output_file, int dim, int nrms, int *rmsindex, double **results, int pixels, int ncustomvalues, char **customnames, int nprintf, int *printfindex, ang_info *angles) {
    printf("\n\n  Writing Results in Output File...\n");
    // Header
    write_2D_control_params_result_matrix_in_file(output_file, NULL, NULL, 1);
    write_attractor_results_matrix_in_file(output_file, NULL, NULL, 1);
    write_min_max_result_matrix_in_file(output_file, dim, NULL, NULL, 1);
    write_min_max_angles_result_matrix_in_file(output_file, angles, NULL, NULL, 1);
    write_rms_result_matrix_in_file(output_file, nrms, rmsindex, NULL, NULL, 1);
    write_customcalc_result_matrix_in_file(output_file, ncustomvalues, nprintf, printfindex, customnames, NULL, NULL, 1);
    fprintf(output_file, "\n");
    // Write Results
    for (int i = 0; i < pixels; i++) {
        int col_offset = 0;
        write_2D_control_params_result_matrix_in_file(output_file, &col_offset, results[i], 2);
        write_attractor_results_matrix_in_file(output_file, &col_offset, results[i], 2);    
        write_min_max_result_matrix_in_file(output_file, dim, &col_offset, results[i], 2);
        write_min_max_angles_result_matrix_in_file(output_file, angles, &col_offset, results[i], 2);
        write_rms_result_matrix_in_file(output_file, nrms, rmsindex, &col_offset, results[i], 2);
        write_customcalc_result_matrix_in_file(output_file, ncustomvalues, nprintf, printfindex, NULL, &col_offset, results[i], 2);
        fprintf(output_file, "\n");
    }
}

void HOS_p_write_fdyndiag_results(FILE *output_file, int dim, int nrms, int *rmsindex, double **results, int pixels, int ncustomvalues, char **customnames, int nprintf, int *printfindex, ang_info *angles) {
    printf("\n\n  Writing Results in Output File...\n");
    // Header
    write_2D_control_params_result_matrix_in_file(output_file, NULL, NULL, 1);
    write_attractor_results_matrix_in_file(output_file, NULL, NULL, 1);
    write_LE_result_matrix_in_file(output_file, dim, NULL, NULL, 1);
    write_min_max_result_matrix_in_file(output_file, dim, NULL, NULL, 1);
    write_min_max_angles_result_matrix_in_file(output_file, angles, NULL, NULL, 1);
    write_rms_result_matrix_in_file(output_file, nrms, rmsindex, NULL, NULL, 1);
    write_customcalc_result_matrix_in_file(output_file, ncustomvalues, nprintf, printfindex, customnames, NULL, NULL, 1);
    fprintf(output_file, "\n");
    // Write Results
    for (int i = 0; i < pixels; i++) {
        int col_offset = 0;
        write_2D_control_params_result_matrix_in_file(output_file, &col_offset, results[i], 2);
        write_attractor_results_matrix_in_file(output_file, &col_offset, results[i], 2);
        write_LE_result_matrix_in_file(output_file, dim, &col_offset, results[i], 2);
        write_min_max_result_matrix_in_file(output_file, dim, &col_offset, results[i], 2);
        write_min_max_angles_result_matrix_in_file(output_file, angles, &col_offset, results[i], 2);
        write_rms_result_matrix_in_file(output_file, nrms, rmsindex, &col_offset, results[i], 2);
        write_customcalc_result_matrix_in_file(output_file, ncustomvalues, nprintf, printfindex, NULL, &col_offset, results[i], 2);
        fprintf(output_file, "\n");
    }
}

