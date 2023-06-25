#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "libs/msg.h"
#include "libs/iofiles.h"
#include "libs/basic.h"

#define BUFSIZE 512

/* Struct Definition and Related Functions */

typedef struct {
    char name[BUFSIZE];
    void *value;
    bool read;
} input;

void init_input_struct(size_t length, input params[length], char *inputnames[]) {
    // Check pointer args
    if (inputnames == NULL) {
        print_debug(" Passing NULL pointer to init_input_struct() function\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }    
    // Loop through all struct array spaces and initialize
    for (int i = 0; i < length; i++) {
        strcpy(params[i].name, inputnames[i]);
        params[i].value = NULL;
        params[i].read = false;
    }
}

/* Safety Check Functions */

void check_for_empty_parameter(char *key, char *svalue) {
    if (strcmp(svalue, "\n") == 0) {
        print_error("INPUT ERROR: parameter '%s' in the input file is empty.\n", key);
        print_error("Add parameter(s) to the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    else {
        return;
    }
}

void check_if_str_is_valid(char *key, char *svalue, char *type, bool only_positive) {
    bool isnumber = check_if_string_is_number(svalue, type, only_positive);
    if (isnumber == false) {
        if (only_positive == true) {
            print_error("INPUT ERROR: '%s' parameter declared in the input file is invalid. It must be a positive '%s' number.\n", key, type);
        }
        else {
            print_error("INPUT ERROR: '%s' parameter declared in the input file is invalid. It must be a '%s' number.\n", key, type);
        }
        print_error("Add the parameter properly before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    else {
        return;
    }
}

void check_for_missing_parameters(input *params, int nparams) {
    // Variables to determine if there is missing parameters in the input file
    bool missing_par = false;
    // Sweep the values of structures to check if there is -1 (not read) in any of the values
    for (int i = 0; i < nparams; i++) {
        if (params[i].read == false) {
            missing_par = true;
            print_error("INPUT ERROR: Cannot find parameter '%s' in the input file.\n", params[i].name);
        }
    }
    // If there is any missing parameter, exit program
    if (missing_par == true) {
        print_error("Add parameter(s) to the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    else {
        return;
    }
}

void check_if_list_exists(char *string, char *defaultstring, char *listname) {
    if (strcmp(string, defaultstring) == 0) { 
        print_error("INPUT ERROR: Cannot find parameter '%s' in the input file.\n", listname);
        print_error("Add parameter's list to the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

void check_for_list_overflow(int n, int listlength, char *listname, char *lengthname) {
    if (n > listlength) {
        print_error("INPUT ERROR: List of '%s' parameters is bigger than '%s' declared in the input file.\n", listname, lengthname);
        print_error("Check parameter(s) in the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

void check_for_index_overflow(input *params, int nparams, char *buffer) {
    // Get the parameter prefix of the array by scanning all the characters 
    // of the first name of struct params until it gets to the character "["
    char param_prefix[strlen(params[0].name) + 1];
    sscanf(params[0].name, "%[^[]", param_prefix);
    // Declare variable to check if there is any index greater than nparams in the input file
    int index = 0;
    if (sscanf(buffer, "%s[^=] = %*f%*[^\n]", param_prefix) == 1) {
        if (sscanf(param_prefix, "%*[^0123456789]%d", &index) == 1) {
            if (index > nparams - 1) {
            print_error("INPUT ERROR: Parameter '%s' is greater than '%d' as declared in the input file.\n", param_prefix, nparams);
            print_error("Check parameter(s) in the input file before running the program again.\n");
            print_exit_prog();
            exit(EXIT_FAILURE);
            }
        }           
    }
}

/* Intermediary Steps for Calling Parameter Reading Functions */

int get_vector_dimension(int n, input *par, char *parname) {
    int value = 0;
    for (int i = 0; i < n; i++) {
        if (strcmp(par[i].name, parname) == 0) {  
            value = *(int*)par[i].value;
        }
    }
    return value;
}

char **define_list_element_names(int nelements, char *basename) {
    // Count dim digits to form "x[dim]" string    
    int digits = count_int_digits(nelements);
    // Allocate memory for a string array that will hold the name "x[0]", "x[1]", "x[2]", ...
    char **names = malloc_string_array(nelements, strlen(basename) + digits + 2 + 1);
    // Copy formatted names to the string array
    for (int i = 0; i < nelements; i++) {
        sprintf(names[i], "%s[%d]", basename, i);
    }
    return names;
}

/* Parameter Reading Functions */

void read_parameter(input *params, char *key, char *svalue, char *type, bool only_positive) {
    // If key is equal to attribute name, then..
    if (strncmp(key, (*params).name, BUFSIZE) == 0) {
        //printf("svalue = %s\n", svalue);
        check_for_empty_parameter(key, svalue);
        // Define the value of the parameter as the stringvalue of the string read
        if (strcmp(type, "int") == 0) {
            check_if_str_is_valid(key, svalue, "int", only_positive);      // Check if string is an int number
            int *tmp_value = malloc(sizeof tmp_value);
            (*tmp_value) = atoi(svalue);
            params->value = tmp_value;
            //printf("key: %s | value: %d\n", key, *(int*)(params->value));
        }
        else if (strcmp(type, "float") == 0) {
            check_if_str_is_valid(key, svalue, "float", only_positive);    // Check if string is a float number
            float *tmp_value = malloc(sizeof tmp_value);
            (*tmp_value) = atof(svalue); 
            params->value = tmp_value;
            //printf("key: %s | value: %f\n", key, *(float*)(params->value));
        }
        else if (strcmp(type, "double") == 0) {
            check_if_str_is_valid(key, svalue, "double", only_positive);   // Check if string is a double number
            double *tmp_value = malloc(sizeof tmp_value);
            (*tmp_value) = atof(svalue); 
            params->value = tmp_value;
            //printf("key: %s | value: %f\n", key, *(double*)(params->value));
        }
        else if ((strcmp(type, "char") == 0) || (strcmp(type, "string") == 0)) {
            size_t len = strlen(svalue);
            params->value = malloc(len + 1);
            strcpy((*params).value, svalue);
            //printf("key: %s | value: %s\n", key, (char*)(params->value));
        }    
        else {
            print_debug("type = '%s' chosen in read_parameter() function. Please, choose one of the valid types: 'int', 'float', 'double', 'char' or 'string'\n", type);
            print_exit_prog();
            exit(EXIT_FAILURE);
        }
        // Flag the value as read (read = true, not read = false)
        params->read = true;

    }

    else {
        return;
    }
}

void read_and_check_parameters(FILE *file, int nparams, input *params, char *param_type, bool only_positive) {
    // Safety check
    if (file == NULL) {
        print_error("File to check parameters could not be openeed.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Set the position to the beggining of the input file
    rewind(file);
    // Allocate memory for buffer
    size_t bufsize = get_size_of_longest_line(file);
    char *buffer = malloc(bufsize);
    rewind(file);
    // Declare variable to hold the value of program parameters read
    int allread;
    while(fgets(buffer, bufsize, file)) {
        // Split line into the buffer with the following separators
        char *key = strtok(buffer, " ,:;={}");
        char *svalue = strtok(NULL, " ,:;={}");
        // Reset allread variable
        allread = 0;
        // Sweep accross all parameters
        for (int i = 0; i < nparams; i++) {
            read_parameter(&params[i], key, svalue, param_type, only_positive);
            // Add a value to allread if the parameter is already read
            if (params[i].read == true) {
                allread++;
            }
        }
        // If allread is equal to the number of parameters, break the loop
        if (allread == nparams) {
            break;
        }
    }
    // Check if all program input parameters were inserted
    check_for_missing_parameters(params, nparams);
    // Free memory
    free(buffer);
}

char *read_and_check_list(FILE *file, char *keyword, char *listname) {
    // Safety check
    if (file == NULL) {
        print_error("File to check parameter list could not be openeed.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Return the position to the start of the file
    rewind(file);
    // Allocate memory for buffer
    size_t bufsize = get_size_of_longest_line(file);
    char *buffer = malloc(bufsize);
    rewind(file);
    // Allocate memory for string
    size_t strsize = get_file_size(file);
    rewind(file);
    // Declare the default string
    char *defaultstring = "empty";
    // Create space in string if strsize if too small
    if (strsize < strlen(defaultstring) + 1) {
        strsize = strlen(defaultstring) + 1;
    }
    // Allocate memory for string
    char *string = malloc(strsize);
    // fill string with the defaultstring ("empty")
    strcpy(string, defaultstring);
    // Declare a variable to store characters from the input file
    char ch;
    // Declare a variable to store the keyword length
    int keywordlength = strlen(keyword);
    // Sweep the input file line by line until it reacher the keyword
    while(fgets(buffer, bufsize, file)) {
        // Check when the buffer starts with keyword
        if (memcmp(buffer, keyword, keywordlength) == 0) {
            // Reset string
            strcpy(string, "");
            // Concatenate the buffer to string
            strcat(string, buffer);
            // Concatenate the remaining of the string until find the enclosing character '}'
            if ((strrchr(string, '}') == NULL)) {  // Check if the 
                while((ch = fgetc(file)) != '}') {
                    strncat(string, &ch, 1);            
                }
                break;
            }
        }
    }

    // Check if string is still the same. If it is, accuse error
    check_if_list_exists(string, defaultstring, listname);
    return string;
}

void read_and_check_list_params(char *list_string, int listsize, input *paramlist, char *listname, char* sizename, char *param_type, bool only_positive) {
    // Declare svalue and get the first token "parameter name" to be discarted
    char *svalue = strtok(list_string, " ,;:={}\n\r\t");
    // Get IC values
    int k = 0;
    // Sweep list_string to get values
    while (svalue != NULL) {
        check_for_list_overflow(k, listsize, listname, sizename);
        // Get next value
        svalue = strtok(NULL, " ,;:={}\n\r\t");
        // Prevent Segmentation fault error if user forget to declare one or more parameters
        if (svalue != NULL) {
            read_parameter(&paramlist[k], paramlist[k].name, svalue, param_type, only_positive);
        }    
        k = k + 1;
    }
    // Check if IC parameter is missing in the input file
    check_for_missing_parameters(paramlist, listsize);
}

void read_and_check_array_parameters(FILE *file, int nparams, input *params, char *param_type, bool only_positive) {
    // Safety check
    if (file == NULL) {
        print_error("File to check parameters could not be openeed.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Set the position to the beggining of the input file
    rewind(file);
    // Allocate memory for buffer
    size_t bufsize = get_size_of_longest_line(file);
    char *buffer = malloc(bufsize);
    rewind(file);
    
    // Read each line of the input file
    while(fgets(buffer, bufsize, file) != NULL) {
        // Split line into the buffer with the following separators
        char *key = strtok(buffer, " ,:;={}");
        char *svalue = strtok(NULL, " ,:;={}");
        // Check if there is any parameter with index greater than nparams in the input file
        check_for_index_overflow(params, nparams, buffer);
        // Sweep accross all parameters and read parameters
        for (int i = 0; i < nparams; i++) {
            read_parameter(&params[i], key, svalue, param_type, only_positive);
        }
    }
    // Check if all program input parameters were inserted
    check_for_missing_parameters(params, nparams);
    // Free memory
    free(buffer);
}

int main(void) {
    // Load File
    char *input_filename = get_input_filename();
    FILE *file = open_file(input_filename, "r", true);
    /* ====================================================================== */
    // HANDLE PROGRAM PARAMETERS
    /* ====================================================================== */
    // Declare parameter names
    char *progparnames[] = {"dim", "npar", "np", "ndiv", "trans", "maxper", 
                            "bifmode", "bifpar"}; 
    // Get the number of parameters in progparnames
    size_t nprogpar = sizeof(progparnames) / sizeof(progparnames[0]);
    // Declare Structure and Initialize it
    input progpar[nprogpar]; 
    init_input_struct(nprogpar, progpar, progparnames);
    // Read prog parameters and check if all were inserted in the input file
    read_and_check_parameters(file, nprogpar, progpar, "int", true);
    /* ====================================================================== */
    // HANDLE INITIAL CONDITIONS
    /* ====================================================================== */
    // Seartch for the value that holds the length of the IC vector
    int dim = get_vector_dimension(nprogpar, progpar, "dim");
    // Declare Struct and variables to store IC parameters and t0
    input ICpar[dim];   // Initial conditions
    input tpar;         // Initial Time
    // Declare parameter names
    char *tname = "t0";
    //char **ICnames = define_IC_names(dim);
    char **ICnames = define_list_element_names(dim, "IC");
    // Initialize structures
    init_input_struct(dim, ICpar, ICnames);
    init_input_struct(1, &tpar, &tname);
    // Read file to store IC parameter list and check if it is inserted in the input file
    char *ICstring = read_and_check_list(file, "IC", "IC = { IC[0], IC[1], ... }");
    read_and_check_list_params(ICstring, dim, ICpar, "IC", "dim (dimension of the system)", "double", true);  
    // Read file to store time parameters and check if it is inserted in the input file
    read_and_check_parameters(file, 1, &tpar, "double", true);
    /* ====================================================================== */
    // HANDLE SYSTEM PARAMETERS
    /* ====================================================================== */
    // Search for the value that holds the number of parameters
    int npar = get_vector_dimension(nprogpar, progpar, "npar");
    // Declare struct to store system parameters
    input syspar[npar];
    // Declare parameter names
    char **sysparnames = define_list_element_names(npar, "p");
    // Initialize structures
    init_input_struct(npar, syspar, sysparnames);
    // Read file to store system parameters and check if it is inserted in the input file
    read_and_check_array_parameters(file, npar, syspar, "double", false);
    /* ====================================================================== */
    // PRINT TO CHECK
    /* ====================================================================== */
    print_warning("\nPROGRAM PARAMETERS:\n");
    for (int i = 0; i < nprogpar; i++) {
        printf("%s = %d\n", progpar[i].name, *(int*)progpar[i].value);
    }
    print_warning("\nINITIAL CONDITIONS:\n");
    for (int i = 0; i < dim; i++) {
        printf("%s = %g\n", ICpar[i].name, *(double*)ICpar[i].value);
    }
    printf("%s = %g\n", tpar.name, *(double*)tpar.value);
    print_warning("\nSYSTEM PARAMETERS:\n");
    for (int i = 0; i < npar; i++) {
        printf("%s = %g\n", syspar[i].name, *(double*)syspar[i].value);
    }


    // NEED CHECKING FOR DUPLICATES!!!


    // free memory
    free(input_filename);
    free_2D_mem((void**)ICnames, dim);
    free_2D_mem((void**)sysparnames, npar);
}