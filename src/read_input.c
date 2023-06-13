#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "libs/msg.h"
#include "libs/iofiles.h"
#include "libs/basic.h"

#define BUFSIZE 512

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

void read_parameter(input *params, char *key, char *svalue, char *type) {
    // If key is equal to attribute name, then..
    if (strncmp(key, (*params).name, BUFSIZE) == 0) {
        printf("svalue = %s\n", svalue);
        check_for_empty_parameter(key, svalue);
        // Define the value of the parameter as the stringvalue of the string read
        if (strcmp(type, "int") == 0) {
            int *tmp_value = malloc(sizeof tmp_value);
            (*tmp_value) = atoi(svalue);
            params->value = tmp_value;
            //printf("key: %s | value: %d\n", key, *(int*)(params->value));
        }
        else if ((strcmp(type, "float") == 0) || (strcmp(type, "double") == 0)) {
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

void read_and_check_prog_parameters(FILE* file, int nprogpar, input *progpar) {
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
        for (int i = 0; i < nprogpar; i++) {
            read_parameter(&progpar[i], key, svalue, "int");
            // Add a value to allread if the parameter is already read
            if (progpar[i].read == true) {
                allread++;
            }
            // Check if the string is a number
            bool isnumber = check_if_string_is_number(svalue, "int");
            if (isnumber == false) {
                print_error("'%s' parameter declared in the input file is invalid. It must be a integer number.\n", key);
                print_error("Add the parameter properly before running the program again.\n");
                print_exit_prog();
                exit(EXIT_FAILURE);
            }
        }
        // If allread is equal to the number of parameters, break the loop
        if (allread == nprogpar) {
            break;
        }
    }
    // Check if all program input parameters were inserted
    check_for_missing_parameters(progpar, nprogpar);
}

int main(void) {
    // Load File
    char *input_filename = get_input_filename();
    FILE *file = open_file(input_filename, "r", true);
    /* HANDLE PROGRAM PARAMETERS */
    // Declare variable names
    char *progparnames[] = {"dim", "npar", "np", "ndiv", "trans", "maxper", 
                            "bifmode", "bifpar"};
    // Get the number of parameters in progparnames
    size_t nprogpar = sizeof(progparnames) / sizeof(progparnames[0]);
    // Declare Structure and Initialize it
    input progpar[nprogpar]; 
    init_input_struct(nprogpar, progpar, progparnames);
    // Read prog parameters and check if all were inserted in the input file
    read_and_check_prog_parameters(file, nprogpar, progpar);

    for (int j = 0; j < nprogpar; j++) {
        printf("progpar[%d] = %s = %d\n", j, progpar[j].name, *(int*)progpar[j].value);
    }

    // PROBLEM BY CHECKING IF A STRING IS A NUMBER BECAUSE SOME STRINGS HAVE "\n" characters in bool check_if_string_is_number(const char *str, const char *type) in basic.c


    // free memory
    free(input_filename);
}