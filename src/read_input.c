#define _CRT_SECURE_NO_WARNINGS
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

void check_if_str_is_valid(char *key, char *svalue, char *type) {
    bool isnumber = check_if_string_is_number(svalue, type, true);
    if ((isnumber == false)) {
        print_error("INPUT ERROR: '%s' parameter declared in the input file is invalid. It must be a positive '%s' number.\n", key, type);
        print_error("Add the parameter properly before running the program again.\n");
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
        //printf("svalue = %s\n", svalue);
        check_for_empty_parameter(key, svalue);
        // Define the value of the parameter as the stringvalue of the string read
        if (strcmp(type, "int") == 0) {
            check_if_str_is_valid(key, svalue, "int");      // Check if string is an int number
            int *tmp_value = malloc(sizeof tmp_value);
            (*tmp_value) = atoi(svalue);
            params->value = tmp_value;
            //printf("key: %s | value: %d\n", key, *(int*)(params->value));
        }
        else if (strcmp(type, "float") == 0) {
            check_if_str_is_valid(key, svalue, "float");    // Check if string is a float number
            float *tmp_value = malloc(sizeof tmp_value);
            (*tmp_value) = atof(svalue); 
            params->value = tmp_value;
            //printf("key: %s | value: %f\n", key, *(float*)(params->value));
        }
        else if (strcmp(type, "double") == 0) {
            check_if_str_is_valid(key, svalue, "double");   // Check if string is a double number
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
        }
        // If allread is equal to the number of parameters, break the loop
        if (allread == nprogpar) {
            break;
        }
    }
    // Check if all program input parameters were inserted
    check_for_missing_parameters(progpar, nprogpar);
}

int search_for_vector_dimension(int n, input *par, char *parname) {
    int value = 0;
    for (int i = 0; i < n; i++) {
        if (strcmp(par[i].name, parname) == 0) {  
            value = *(int*)par[i].value;
        }
    }
    return value;
}

void check_for_empty_list(char *string, char *defaultstring, char *listname) {
    if (strcmp(string, defaultstring) == 0) { 
        print_error("INPUT ERROR: Cannot find parameter '%s' in the input file.\n", listname);
        print_error("Add parameter(s) to the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

char *read_and_check_list(FILE *file, char *keyword, char *listname) {
    // Return the position to the start of the file
    rewind(file);
    // Allocate memory for buffer
    size_t bufsize = get_size_of_longest_line(file);
    char *buffer = malloc(bufsize);
    rewind(file);
    // Allocate memory for string
    size_t strsize = get_file_size(file);
    rewind(file);
    // Declare string and fill it with the tag empty
    char *string = malloc(strsize);
    char *defaultstring = "empty";
    strcat(string, defaultstring);
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
    check_for_empty_list(string, defaultstring, listname);
    return string;
}



int main(void) {
    // Load File
    char *input_filename = get_input_filename();
    FILE *file = open_file(input_filename, "r", true);
    /* ====================================================================== */
    // HANDLE PROGRAM PARAMETERS
    /* ====================================================================== */
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
    /* ====================================================================== */
    // HANDLE INITIAL CONDITIONS
    /* ====================================================================== */
    // Seartch for the value that holds the length of the IC vector
    int dim = search_for_vector_dimension(nprogpar, progpar, "dim");
    // Declare Struct and variables to store IC parameters and t0
    input ICpar[dim];   // Initial conditions
    input tpar;         // Initial Time
    char *tname = "t0";
    // Read file to store IC parameter list and check if it is inserted in the input file
    char *ICstring = read_and_check_list(file, "IC", "IC = {x[0], x[1], ...}");

    // PROBLEM IN CHECKING AN EMPTY LIST in check_for_empty_list(string, defaultstring, listname)


    printf("ICstring: %s\n", ICstring);
    
    
    for (int j = 0; j < nprogpar; j++) {
        printf("progpar[%d] = %s = %d\n", j, progpar[j].name, *(int*)progpar[j].value);
    }
    

    // free memory
    free(input_filename);
}