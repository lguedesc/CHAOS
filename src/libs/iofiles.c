#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#ifdef _WIN32
    #include <direct.h>
    const char SEP = '\\';
#else
    const char SEP = '/';
#endif    

// Define the maximum number of repeatable files that the program can create in the directory
#define MAX_FILE_NUM 1000
// Define function to print name of the variable in a string
#define getName(var)  #var

struct stat info;

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
    int n = 200;                                // n*sizeof(char) = max size of path
    // Declare variables
    const size_t len = strlen(pathname);        // Size of the pathname
    char *path = malloc(n * sizeof(char));    // Allocate temporary memory to copy pathname
    char *ch;                                   // Pointer to each character in the string
    // Check if len is too long to be stored into path
    if (len > n * sizeof(char) - 1) {
        errno = ENAMETOOLONG;                   // Return error if its too long
        perror(pathname);                       // Print error
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
        printf("Cannot find output directory\nCreating output directory...\n");
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

char *get_input_filename(void) {
    char *filename = malloc(200 * sizeof(*filename));
    printf("  Enter Input Filename: ");
    scanf("%s", filename);
    return filename;
    /* The user is responsible to free (filename) after the function call */
}

// Write Results

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

void write_bifurc_results(FILE *output_file, int dim, double varpar, double *x, int mode) {
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
    else {
        printf("Failed to write results in output file using mode (%d)...\n", mode);
        return;
    }
}

void write_fbifurc_results(FILE *output_file, int dim, int np, int trans, double varpar, double *x, int attractor, double **poinc, int diffattrac, int mode) {
    // Check the mode of the function
    // Header
    if (mode == 0) {
        fprintf(output_file, "%s ", "Cpar");
        for (int i = 0; i < dim; i++) {
            fprintf(output_file, "x[%i] ", i);
        }
        fprintf(output_file, "Attractor DiffAttrac\n");
    } 
    // Write Results
    else if (mode == 1) {
        for(int q = 0; q < np - trans; q ++) {
            fprintf(output_file, "%.10lf ", varpar);
            for (int i = 0; i < dim; i++) {
                fprintf(output_file, "%.10lf ", poinc[q][i]);
            }
            fprintf(output_file, "%i %i\n", attractor, diffattrac);
        }
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

void p_write_dyndiag_results(FILE *output_file, int dim, double **results, int pixels) {
    // Header
    fprintf(output_file, "%s %s %s %s ", "CparY", "CparX", "Attractor", "DiffAttrac");
    for (int i = 0; i < dim; i++) {
        fprintf(output_file, "LE[%d] ", i);
    }
    fprintf(output_file, "\n");
    // Write Results
    for (int i = 0; i < pixels; i++) {
        fprintf(output_file, "%.10lf %.10lf %d %d ", results[i][0], results[i][1], (int)results[i][2], (int)results[i][3]);
        for (int j = 0; j < dim; j++) {
            fprintf(output_file, "%.10lf ", results[i][j+4]);
        }
        fprintf(output_file, "\n");
    }
}

void p_write_epbasin_results(FILE *output_file, double **results, int pixels, int dim) {
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

