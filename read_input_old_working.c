#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "src/libs/msg.h"
#include "src/libs/iofiles.h"
#include "src/libs/basic.h"

#define BUFSIZE 512

/* Struct Definition and Related Functions */

typedef struct {
    char name[BUFSIZE];
    void *value;
    bool read;
} input;

static void init_input_struct(size_t length, input params[length], char *inputnames[]) {
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

void free_input_struct(size_t length, input *params) {
    for (int i = 0; i < length; i++) {
        free(params[i].value);
    }
    free(params);
}

/* Safety Check Functions */

static void check_for_empty_parameter(char *key, char *svalue) {
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

static void check_if_str_is_valid(char *key, char *svalue, char *type, bool only_positive) {
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

static void check_for_missing_parameters(input *params, int nparams) {
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

static void check_if_list_exists(char *string, char *defaultstring, char *listname) {
    if (strcmp(string, defaultstring) == 0) { 
        print_error("INPUT ERROR: Cannot find parameter '%s' in the input file.\n", listname);
        print_error("Add parameter's list to the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

static void check_for_list_overflow(int n, int listlength, char *listname, char *lengthname) {
    if (n > listlength) {
        print_error("INPUT ERROR: List of '%s' parameters is bigger than '%s' declared in the input file.\n", listname, lengthname);
        print_error("Check parameter(s) in the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

static void check_for_index_overflow(input *params, int nparams, char *buffer) {
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

static void check_for_duplicate_parameter(FILE *file, input *params, int nparams, char *buffer, size_t bufsize) {
    // Safety check
    file_safety_check(file);
    // Declare variable to check how many times a name was found
    int found;
    // Declare variable to check if programs needs to be closed
    bool close = false;
    // Loop through each name in the struct params
    for (int i = 0; i < nparams; i++) {
        // Reset value of found
        found = 0;
        // Go to the beggining of the file
        rewind(file);
        // Read line by line searching for keywords
        while(fgets(buffer, bufsize, file) != NULL) {
            // Split line into the buffer with the following separators
            // to get the keyword
            char *key = strtok(buffer, " ,:;={}");
            // Check if key is equal to param_names[i]
            if (strcmp(key, params[i].name) == 0) {
                // Increment found if it is equal
                found++;
            }
        }
        // Check if name was found more than once
        if (found > 1) {
            close = true;
            print_error("%d instances of parameter '%s' found within the input file.\n", found, params[i].name);
        }
    }
    // Check if the there is duplicated to kill execution
    if (close == true) {
        print_error("Check parameters in the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    } else {
        return;
    }
}  

static void check_for_duplicate_list(FILE *file, char* keyword, char *buffer, size_t bufsize) {
    // Safety check
    file_safety_check(file);
    // Declare variable to check how many times the listname was found
    int found = 0;
    // Go to the beggining of the file
    rewind(file);
    // Read line by line searching for keyword: listname
    while (fgets(buffer, bufsize, file) != NULL) {
        // Split line into the buffer with the following separators
        // to get the keyword
        char *key = strtok(buffer, " ,:;={}");
        // Check if key is equal to listname
        if (strcmp(key, keyword) == 0) {
            found++;
        }
    }
    if (found > 1) {
        print_error("%d instances of list '%s' found within the input file.\n", found, keyword);
        print_error("Check parameters in the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    } else {
        return;
    }
}

static bool check_optional_parameter(FILE *file, char *param_name) {
    bool declared = false;
    // Safety check
    file_safety_check(file);
    // Set the position to the beggining of the input file
    rewind(file);
    // Allocate memory for buffer
    size_t bufsize = get_size_of_longest_line(file);;
    char *buffer = malloc(bufsize);
    // Set the position to the beggining of the input file again
    rewind(file);
    // Read the file line by line
    while(fgets(buffer, bufsize, file) != NULL) {
        // Split line into the buffer with the following separators
        char *key = strtok(buffer, " ,:;={}");
        if (strcmp(key, param_name) == 0) {
            declared = true;
        }
    }
    // Free memory
    free(buffer);

    return declared;
}

static void check_out_of_range_list_param(int limit, input *param, int nparams) {
    // Variable to check if there is a invalid value
    bool invalid = true;
    for (int i = 0; i < nparams; i++) {
        if (*(int*)param[i].value > limit - 1) {
            print_error("Parameter '%s' is out of range. The list must contain parameters from 0 to %d.\n", param[i].name, limit-1);
            print_error("Check parameters in the input file before running the program again.\n");
            print_exit_prog();
            exit(EXIT_FAILURE);
        }
    }
}

/* Intermediary Steps for Calling Parameter Reading Functions */
/*
static int get_vector_dimension(int n, input *par, char *parname) {
    int value = 0;
    for (int i = 0; i < n; i++) {
        if (strcmp(par[i].name, parname) == 0) {  
            value = *(int*)par[i].value;
        }
    }
    return value;
}*/

static int find_param_value_by_name(input *par, char *parname) {
    int value = 0;
    int i = 0;
    while (par[i].name[0] != '\0') {
        if (strcmp(par[i].name, parname) == 0) {
            value = *(int*)par[i].value;
            break;
        }
        i++;
    }
    return value;
}

static char **define_list_element_names(int nelements, char *basename) {
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

static void read_parameter(input *params, char *key, char *svalue, char *type, bool only_positive) {
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

static void read_and_check_parameters(FILE *file, int nparams, input *params, char *param_type, bool only_positive) {
    // Safety check
    file_safety_check(file);
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

static char *read_and_check_list(FILE *file, char *keyword, char *listname) {
    // Safety check
    file_safety_check(file);
    // Return the position to the start of the file
    rewind(file);
    // Allocate memory for buffer
    size_t bufsize = get_size_of_longest_line(file);
    char *buffer = malloc(bufsize);
    rewind(file);
    // Check for list duplicate
    check_for_duplicate_list(file, keyword, buffer, bufsize);
    // Set the position to the beggining of the input file again
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

static void read_and_check_list_params(char *list_string, int listsize, input *paramlist, char *listname, char* sizename, char *param_type, bool only_positive) {    
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

static void read_and_check_array_parameters(FILE *file, int nparams, input *params, char *param_type, bool only_positive) {
    // Safety check
    file_safety_check(file);
    // Set the position to the beggining of the input file
    rewind(file);
    // Allocate memory for buffer
    size_t bufsize = get_size_of_longest_line(file);
    char *buffer = malloc(bufsize);
    check_for_duplicate_parameter(file, params, nparams, buffer, bufsize);
    // Set the position to the beggining of the input file again
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

/* Functions to read specific block of parameters */
static void handle_program_params(FILE *file, size_t nprogpar, input *progpar, char **progparnames) {
    // Safety check
    file_safety_check(file);
    // Initialize struct
    init_input_struct(nprogpar, progpar, progparnames);
    // Read prog parameters and check if all were inserted in the input file
    read_and_check_parameters(file, nprogpar, progpar, "int", true);
}

static void handle_initial_conditions(FILE *file, int dim, input *ICpar, input *tpar) {
    // Declare parameter names
    char *tname = "t0";
    char **ICnames = define_list_element_names(dim, "IC");
    // Initialize structures
    init_input_struct(dim, ICpar, ICnames);
    init_input_struct(1, tpar, &tname);
    // Read file to store IC parameter list and check if it is inserted in the input file
    char *ICstring = read_and_check_list(file, "IC", "IC = { IC[0], IC[1], ... }");
    read_and_check_list_params(ICstring, dim, ICpar, "IC", "dim (dimension of the system)", "double", true);  
    // Read file to store time parameters and check if it is inserted in the input file
    read_and_check_parameters(file, 1, tpar, "double", true);
    // Free memory
    free_2D_mem((void**)ICnames, dim);
}

static void handle_system_params(FILE *file, int npar, input *syspar) {
    // Safety check
    file_safety_check(file);
    // Declare parameter names
    char **sysparnames = define_list_element_names(npar, "p");
    // Initialize structures
    init_input_struct(npar, syspar, sysparnames);
    // Read file to store system parameters and check if it is inserted in the input file
    read_and_check_array_parameters(file, npar, syspar, "double", false);
    // Free memory
    free_2D_mem((void**)sysparnames, npar);
}

static void handle_bifurc_params(FILE *file, int nbifpar, int nbiflims, input *bifpar, input *biflimits, char **bifnames) {
    // Safety check
    file_safety_check(file);
    // Declare parameter names
    char **biflimnames = define_list_element_names(nbiflims, "biflimits");
    // Initialize structs
    init_input_struct(nbifpar, bifpar, bifnames);
    init_input_struct(nbiflims, biflimits, biflimnames);
    // Read file to store bifurcation limits and check if it is properly defined in the input file
    read_and_check_parameters(file, nbifpar, bifpar, "int", true);
    char *bifstring = read_and_check_list(file, "biflimits", "biflimits = { initial, final, number of points }");
    read_and_check_list_params(bifstring, nbiflims, biflimits, "biflimits", "3", "double", false);
    // Free memory
    free_2D_mem((void**)biflimnames, nbiflims);
}

static void handle_optional_list_params(FILE *file, int list_len, input *list, int max_number_element, char *base_listname, char *list_description, char *size_description) {
    // Safety Check
    file_safety_check(file);
    // Define list name
    char **listname = define_list_element_names(list_len, base_listname);
    // Initialize structure
    init_input_struct(list_len, list, listname);
    // Read file to store values of parameter list and check if it is inserted in the input file
    char *liststring = read_and_check_list(file, base_listname, list_description);
    read_and_check_list_params(liststring, list_len, list, base_listname, size_description, "int", true);
    // Check if the elements in list are invalid (greater than max_number_element)
    check_out_of_range_list_param(max_number_element, list, list_len);
    // Free memory
    free_2D_mem((void**)listname, list_len);
}   

/* Functions to assign values */

void get_program_params(FILE *file, input **progpar, size_t *nprogpar) {
    // Declare parameter names
    char *progparnames[] = {"dim", "npar", "np", "ndiv", "trans", "maxper"}; 
    // Get the number of parameters in progparnames
    (*nprogpar) = sizeof(progparnames) / sizeof(progparnames[0]);
    // Declare Structure and Initialize it
    (*progpar) = malloc((*nprogpar) * sizeof(**progpar)); 
    // Get program params
    handle_program_params(file, (*nprogpar), (*progpar), progparnames);
}

void get_initial_conditions(FILE *file, input *t0par, input **ICpar, int dim) {
    // Declare Struct and variables to store IC parameters and t0
    (*ICpar) = malloc(dim * sizeof(**ICpar));   // Initial conditions
    // Get program params
    handle_initial_conditions(file, dim, (*ICpar), t0par);
}

void get_system_params(FILE *file, input **syspar, int npar) {
    // Declare struct to store system parameters
    (*syspar) = malloc(npar * sizeof(**syspar));
    // Get system's parameters
    handle_system_params(file, npar, (*syspar));
}

void get_bifurc_params(FILE *file, input **bifpar, input **biflimits) {
    // Declare parameter names
    char *bifnames[] = {"bifmode", "bifpar"}; 
    // Get the number of parameters for each input struct
    size_t nbifpar = sizeof(bifnames) / sizeof(bifnames[0]);
    int nbiflims = 3;
    // Declare Structures 
    (*bifpar) = malloc(nbifpar * sizeof(**bifpar));
    (*biflimits) = malloc(nbiflims * sizeof(**biflimits));
    // Get bifurcation parameters
    handle_bifurc_params(file, nbifpar, nbiflims, (*bifpar), (*biflimits), bifnames);
}

void get_rms_params(FILE *file, input *nrms, input **rmslist, int dim) {
    // Declare the name of the control parameter
    char *rmsname = "nrms";
    init_input_struct(1, nrms, &rmsname);
    // Check if nrms is declared in the input file
    bool rmscalc = check_optional_parameter(file, rmsname);
    if (rmscalc == true) {        
        // Give the values to struct nrms
        read_and_check_parameters(file, 1, nrms, "int", true);
        // Get length of the rms list
        int list_len = *(int*)(*nrms).value;
        // Check if nrms is greater than zero to continue
        if (list_len > 0) {
            // Declare struct to store RMS list
            (*rmslist) = malloc(list_len * sizeof(**rmslist));
            // Get RMS parameters         
            //handle_rms_params(file, list_len, rmslist, dim);
            handle_optional_list_params(file, list_len, (*rmslist), dim, "rms", "rms = { variable index 0, variable index 1, ... }", "nrms (number of rms calculations)");
        }
    } else {
        (*nrms).value = 0;
        (*nrms).read = true;
    }
}

void get_ccalc_params(FILE *file, input *nccalc, input *nend_ccalc, input *nbody_ccalc, input **end_ccalc_list, input **body_ccalc_list) {
    // Declare names of params
    char *nccalc_name = "nccalc";
    char *nend_ccalc_name = "nccalc_";
    char *nbody_ccalc_name = "nccalc*";
    // Initialize structs
    init_input_struct(1, nccalc, &nccalc_name);
    init_input_struct(1, nend_ccalc, &nend_ccalc_name);
    init_input_struct(1, nbody_ccalc, &nbody_ccalc_name);
    // Check if ccalc is declared in the input file
    bool ccalc = check_optional_parameter(file, nccalc_name);
    if (ccalc == true) {
        // Give the values to struct nccalc
        read_and_check_parameters(file, 1, nccalc, "int", true);
        // Get max element
        int max_element = *(int*)(*nccalc).value;
        if (max_element > 0) {
            // Check if end ccalc and body ccalc are declared in the input file
            bool end_ccalc = check_optional_parameter(file, nend_ccalc_name);
            bool body_ccalc = check_optional_parameter(file, nbody_ccalc_name);
            if (end_ccalc == true) {
                // Give values to struct nend_ccalc
                read_and_check_parameters(file, 1, nend_ccalc, "int", true);
                // Get length of nend_ccalc list
                int end_list_len = *(int*)(*nend_ccalc).value;
                if (end_list_len > 0) {
                    // Declare struct to store end_ccalc list
                    (*end_ccalc_list) = malloc(end_list_len * sizeof(**end_ccalc_list));
                    // Get end_ccalc parameters
                    handle_optional_list_params(file, end_list_len, (*end_ccalc_list), max_element, "ccalc_", "ccalc_ = { ccalc_[0], ccalc_[1], ... }", "nccalc (number of custom calculations that will be printed)");
                }
            } 
            else {
                // Put default values and mark as read
                (*nend_ccalc).value = 0;
                (*nend_ccalc).read = true;
            }
            if (body_ccalc == true) {
                // Give values to struct nbody_ccalc
                read_and_check_parameters(file, 1, nbody_ccalc, "int", true);
                // Get length of nbody_ccalc list
                int body_list_len = *(int*)(*nbody_ccalc).value;
                if (body_list_len > 0) {
                    // Declare struct to store body_ccalc list
                    (*body_ccalc_list) = malloc(body_list_len * sizeof(**body_ccalc_list));
                    // Get body ccalc parameters
                    handle_optional_list_params(file, body_list_len, (*body_ccalc_list), max_element, "ccalc*", "ccalc* = { ccalc*[0], ccalc*[1], ... }", "nccalc (number of custom calculations that will be printed)");
                }
            }
            else {
                // Put default values and mark as read
                (*nbody_ccalc).value = 0;
                (*nbody_ccalc).read = true;
            }
        }        
    }
    else {
        // Put default values and mark as read
        (*nccalc).value = 0;
        (*nccalc).read = true;
        (*nend_ccalc).value = 0;
        (*nend_ccalc).read = true;
        (*nbody_ccalc).value = 0;
        (*nbody_ccalc).read = true;
    }
}

void print_to_check(int nprogpar, int dim, int npar, input *progpar, input *ICpar, input t0par, input *syspar, input *bifpar, input *biflimits, input nrms,
                    input *rmslist, input nccalc, input nend_ccalc, input *end_ccalc_list, input nbody_ccalc, input *body_ccalc_list) {
    // Print to check
    print_warning("\nPROGRAM PARAMETERS:\n");
    for (int i = 0; i < nprogpar; i++) {
        printf("%s = %d\n", progpar[i].name, *(int*)progpar[i].value);
    }
    print_warning("\nINITIAL CONDITIONS:\n");
    for (int i = 0; i < dim; i++) {
        printf("%s = %g\n", ICpar[i].name, *(double*)ICpar[i].value);
    }
    printf("%s = %g\n", t0par.name, *(double*)t0par.value);
    print_warning("\nSYSTEM PARAMETERS:\n");
    for (int i = 0; i < npar; i++) {
       printf("%s = %g\n", syspar[i].name, *(double*)syspar[i].value);
    }
    print_warning("\nBIFURCATION PARAMETERS:\n");
    for (int i = 0; i < 2; i++) {
        printf("%s = %d\n", bifpar[i].name, *(int*)bifpar[i].value);
    }
    for (int i = 0; i < 3; i++) {
        printf("%s = %g\n", biflimits[i].name, *(double*)biflimits[i].value);
    }
    if (*(int*)nrms.value > 0) {
        print_warning("\nRMS PARAMETERS:\n");
        printf("%s = %d\n", nrms.name, *(int*)nrms.value);
        for (int i = 0; i < *(int*)nrms.value; i++) {
            printf("%s = %d\n", rmslist[i].name, *(int*)rmslist[i].value);
        }
    }
    if (*(int*)nccalc.value > 0) {
        print_warning("\nCUSTOM CALCULATIONS PARAMETERS:\n");
        printf("%s = %d\n", nccalc.name, *(int*)nccalc.value);
        if (*(int*)nend_ccalc.value > 0) {
            printf("%s = %d\n", nend_ccalc.name, *(int*)nend_ccalc.value);
            for (int i = 0; i < *(int*)nend_ccalc.value; i++) {
                printf("%s = %d\n", end_ccalc_list[i].name, *(int*)end_ccalc_list[i].value);
            }
        } 
        if (*(int*)nbody_ccalc.value > 0) {
            printf("%s = %d\n", nbody_ccalc.name, *(int*)nbody_ccalc.value);
            for (int i = 0; i < *(int*)nbody_ccalc.value; i++) {
                printf("%s = %d\n", body_ccalc_list[i].name, *(int*)body_ccalc_list[i].value);
            }
        }
    }
}

int main(void) {
    // Load File
    char *input_filename = get_input_filename();
    FILE *file = open_file(input_filename, "r", true);
    /* ====================================================================== */
    // HANDLE PROGRAM PARAMETERS
    /* ====================================================================== */
    input *progpar = NULL;
    size_t nprogpar = 0;
    get_program_params(file, &progpar, &nprogpar);
    /* ====================================================================== */
    // HANDLE INITIAL CONDITIONS
    /* ====================================================================== */
    // Seartch for the value that holds the length of the IC vector
    int dim = find_param_value_by_name(progpar, "dim");
    input t0par;
    input *ICpar = NULL;
    get_initial_conditions(file, &t0par, &ICpar, dim);
    /* ====================================================================== */
    // HANDLE SYSTEM PARAMETERS
    /* ====================================================================== */
    // Search for the value that holds the number of parameters
    int npar = find_param_value_by_name(progpar, "npar");
    input *syspar = NULL;
    get_system_params(file, &syspar, npar);    
    /* ====================================================================== */
    // HANDLE BIFURCATION PARAMETERS
    /* ====================================================================== */
    input *bifpar = NULL;
    input *biflimits = NULL;
    get_bifurc_params(file, &bifpar, &biflimits);
    /* ====================================================================== */
    // HANDLE RMS PARAMETERS
    /* ====================================================================== */
    // Declare struct to store RMS list size, initialize it and read from file
    input nrms;
    input *rmslist = NULL;
    get_rms_params(file, &nrms, &rmslist, dim);
    /* ====================================================================== */
    // HANDLE CUSTOMCALC PARAMETERS
    /* ====================================================================== */
    input nccalc;
    input nend_ccalc;
    input nbody_ccalc;
    input *end_ccalc_list = NULL; 
    input *body_ccalc_list = NULL;
    get_ccalc_params(file, &nccalc, &nend_ccalc, &nbody_ccalc, &end_ccalc_list, &body_ccalc_list);
    
    // Print to check
    print_to_check(nprogpar, dim, npar, progpar, ICpar, t0par, syspar, bifpar, biflimits, nrms, rmslist, nccalc, nend_ccalc, end_ccalc_list, nbody_ccalc, body_ccalc_list);
    // Free memory
    free(input_filename);
    free_input_struct(nprogpar, progpar);
    free_input_struct(dim, ICpar);
    free_input_struct(npar, syspar);
    free_input_struct(2, bifpar);
    free_input_struct(3, biflimits);
    free_input_struct(*(int*)nrms.value, rmslist);
    free_input_struct(*(int*)nend_ccalc.value, end_ccalc_list);
    free_input_struct(*(int*)nbody_ccalc.value, body_ccalc_list);
}