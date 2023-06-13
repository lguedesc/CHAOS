#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define BUFSIZE 512

typedef struct {
    char name[BUFSIZE];
    double value;
    int read;
} doublepar;

typedef struct {
    char name[BUFSIZE];
    int value;
    int read;
} intpar;

void read_double_parameter(doublepar *attr, char *key, char *svalue) {
    // If key is equal to attribute name
    if (strncmp(key, (*attr).name, BUFSIZE) == 0) {
        // Define the value of the attribue as the stringvalue of the string read
        (*attr).value = atof(svalue);
        // Flag the value as read (read = 1, not read = -1)
        (*attr).read = 1;
        //printf("key: %s | value: %lf\n", key, (*attr).value);
    } 
    else {
        return;
    }
}

void read_int_parameter(intpar *attr, char *key, char *svalue) {
    // If key is equal to attribute name, then..
    if (strncmp(key, (*attr).name, BUFSIZE) == 0) {
        // Define the value of the attribue as the stringvalue of the string read
        (*attr).value = atoi(svalue);
        // Flag the value as read (read = 1, not read = -1)
        (*attr).read = 1;
        //printf("key: %s | value: %d\n", key, (*attr).value);
    } 
    else {
        return;
    }
}

void check_for_double_missing_parameters(doublepar *doubles, int ndoubles) {
    // Variables to determine if there is missing parameters in the input file
    int missing_par = 0;
    // Sweep the values of structures to check if there is -1 (not read) in any of the values
    for (int i = 0; i < ndoubles; i++) {
        if (doubles[i].read == -1) {
            missing_par = 1;
            printf("  INPUT ERROR: Cannot find parameter '%s' in the input file.\n", doubles[i].name);
        }
    }
    // If there is any missing parameter, exit program
    if (missing_par == 1) {
        printf("               Add parameter(s) to the input file before running the program again.\n");
        printf("               Exiting Program...\n");
        exit(EXIT_FAILURE);
    }
    else {
        return;
    }
}

void check_for_int_missing_parameters(intpar *ints, int nints) {
    // Variables to determine if there is missing parameters in the input file
    int missing_par = 0;
    // Sweep the values of structures to check if there is -1 (not read) in any of the values
    for (int i = 0; i < nints; i++) {
        if (ints[i].read == -1) {
            missing_par = 1;
            printf("  INPUT ERROR: Cannot find parameter '%s' in the input file.\n", ints[i].name);
        }
    }
    // If there is any missing parameter, exit program
    if (missing_par == 1) {
        printf("               Add parameter(s) to the input file before running the program again.\n");
        printf("               Exiting Program...\n");
        exit(EXIT_FAILURE);
    }
    else {
        return;
    }
}

void check_for_input_list(char *string, char *defaultstring, char *listname) {
    if (strcmp(string, defaultstring) == 0) { 
        printf("  INPUT ERROR: Cannot find parameter '%s' in the input file.\n", listname);
        printf("               Add parameter(s) to the input file before running the program again.\n");
        printf("               Exiting Program...\n");
        exit(EXIT_FAILURE);
    }
}

int search_for_vector_dimension(int n, intpar *par, char *parname) {
    int value = 0;
    for (int i = 0; i < n; i++) {
        if (strcmp(par[i].name, parname) == 0) {
            value = par[i].value;   
        }
    }
    return value;
}

char *read_and_check_list(FILE *input, char *buffer, char *keyword, char *listname) {
    // Return the position to the start of the file
    rewind(input);
    // Declare string and fill it with the tag empty
    char *string = malloc(BUFSIZE * sizeof string);
    char *defaultstring = "empty";
    strcat(string, defaultstring);
    // Declare a variable to store characters from the input file
    char ch;
    // Declare a variable to store the keyword length
    int keywordlength = strlen(keyword);
    // Sweep the input file line by line until it reacher the keyword
    while(fgets(buffer, BUFSIZE, input)) {
        // Check when the buffer starts with keyword
        if (memcmp(buffer, keyword, keywordlength) == 0) {
            // Reset string
            strcpy(string, "");
            // Concatenate the buffer to string
            strcat(string, buffer);
            // Concatenate the remaining of the string until find the enclosing character '}'
            if ((strrchr(string, '}') == NULL)    ) {  // Check if the 
                while((ch = fgetc(input)) != '}') {
                    strncat(string, &ch, 1);            
                }
                break;
            }
        }
    }
    // Check if string is still the same. If it is, accuse error
    check_for_input_list(string, defaultstring, listname);
    return string;
}

void check_for_list_overflow(int n, int listlength, char *listname, char *lengthname) {
    if (n > listlength) {
        printf("  INPUT ERROR: List of '%s' parameters is bigger than '%s' declared in the input file.\n", listname, lengthname);
        printf("               Check parameter(s) in the input file before running the program again.\n");
        printf("               Exiting Program...\n");
        exit(EXIT_FAILURE);
    }
}

void check_for_invalid_int_input_value(intpar *ints, int nints, int limit, char *limitname, int mode) {
    // Variables to determine if any declared element
    int wrong_par = 0;
    /* Sweep input values to check if there is invalid entries */
    // case 1: 0 to a certain range limit
    if (mode == 0) {
        for (int i = 0; i < nints; i++) {
            if (ints[i].value > limit) {  // nrms > 5
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is bigger than '%s'. It must be smaller or equal.\n", ints[i].name, limitname);
            }
            else if (ints[i].value < 0) { // nrms < 0
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is negative. It must be positive.\n", ints[i].name);
            }
        }
    } 
    // case 2: only positive values
    else if (mode == 1) {  
        for (int i = 0; i < nints; i++) {
            if (ints[i].value < 0) {    // nrms < 0
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is negative. It must be positive.\n", ints[i].name);
            }
        }
    }
    // case 3: only negative values
    else if (mode == 2) {
        for (int i = 0; i < nints; i++) {
            if (ints[i].value > 0) {  // nrms > 0
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is positive. It must be negative.\n", ints[i].name);
            }
        }
    }
    // case 4: range limit negative value to 0
    else if (mode == 3) { 
        for (int i = 0; i < nints; i++) {
            if (ints[i].value < limit) { // nrms < -3
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is smaller than '%s'. It must be bigger or equal.\n", ints[i].name, limitname);
            }
            else if (ints[i].value > 0) {
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is positive. It must be begative.\n", ints[i].name);
            }
        }
    }
    // If there is any invalid value, exit program
    if (wrong_par == 1) {
        printf("               Check the input file before running the program again.\n");
        printf("               Exiting Program...\n");
        exit(EXIT_FAILURE);
    }
    else {
        return;
    }
}

void check_for_invalid_double_input_value(doublepar *ints, int nints, int limit, char *limitname, int mode) {
    // Variables to determine if any declared element
    int wrong_par = 0;
    /* Sweep input values to check if there is invalid entries */
    // case 1: 0 to a certain range limit
    if (mode == 0) {
        for (int i = 0; i < nints; i++) {
            if (ints[i].value > limit) {  // nrms > 5
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is bigger than '%s'. It must be smaller or equal.\n", ints[i].name, limitname);
            }
            else if (ints[i].value < 0) { // nrms < 0
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is negative. It must be positive.\n", ints[i].name);
            }
        }
    } 
    // case 2: only positive values
    else if (mode == 1) {  
        for (int i = 0; i < nints; i++) {
            if (ints[i].value < 0) {    // nrms < 0
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is negative. It must be positive.\n", ints[i].name);
            }
        }
    }
    // case 3: only negative values
    else if (mode == 2) {
        for (int i = 0; i < nints; i++) {
            if (ints[i].value > 0) {  // nrms > 0
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is positive. It must be negative.\n", ints[i].name);
            }
        }
    }
    // case 4: range limit negative value to 0
    else if (mode == 3) { 
        for (int i = 0; i < nints; i++) {
            if (ints[i].value < limit) { // nrms < -3
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is smaller than '%s'. It must be bigger or equal.\n", ints[i].name, limitname);
            }
            else if (ints[i].value > 0) {
                wrong_par = 1;
                printf("  INPUT ERROR: Parameter '%s' is positive. It must be begative.\n", ints[i].name);
            }
        }
    }
    // If there is any invalid value, exit program
    if (wrong_par == 1) {
        printf("               Check the input file before running the program again.\n");
        printf("               Exiting Program...\n");
        exit(EXIT_FAILURE);
    }
    else {
        return;
    }
}

void initialize_doublepar_list_struct(doublepar *attr, int vectorlength, char *varname) {
    char attrnames[20];
    for (int i = 0; i < vectorlength; i++) {
        sprintf(attrnames, "%s[%d]", varname, i);   // Determine name
        strcpy(attr[i].name, attrnames);    // Insert in name .name
        attr[i].value = 0.0;                // Determine initial value
        attr[i].read = -1;                  // Determine read status 
    }
}

void initialize_doublepar_struct(doublepar *attr, int vectorlength, char *varnames[]) {
    if ((attr == NULL) || (varnames == NULL)) {
        printf("  Passing NULL pointer(s) to initialize_intpar_struct() function");
        return;
    }
    for (int i = 0; i < vectorlength; i++) {
        strcpy(attr[i].name, varnames[i]);
        attr[i].value = 0.0; 
        attr[i].read = -1;
    }
}

void initialize_intpar_struct(intpar *attr, int vectorlength, char *varnames[]) {
    if ((attr == NULL) || (varnames == NULL)) {
        printf("  Passing NULL pointer(s) to initialize_intpar_struct() function");
        return;
    }
    for (int i = 0; i < vectorlength; i++) {
        strcpy(attr[i].name, varnames[i]);
        attr[i].value = 0; 
        attr[i].read = -1;
    }
}

void initialize_intpar_list_struct(intpar *attr, int vectorlength, char *varname) {
    char attrnames[20];
    for (int i = 0; i < vectorlength; i++) {
        sprintf(attrnames, "%s[%d]", varname, i);   // Determine name
        strcpy(attr[i].name, attrnames);    // Insert in name .name
        attr[i].value = 0;                  // Determine initial value
        attr[i].read = -1;                  // Determine read status 
    }
}

void read_and_check_prog_parameters(FILE* input, char *buffer, int nprogpar, intpar *progpar) {
    // Set the position to the beggining of the input file
    rewind(input);
    // Declare variable to hold the value of program parameters read
    int allread;
    while(fgets(buffer, BUFSIZE, input)) {
        // Split line into the buffer with the following separators
        char *key = strtok(buffer, " ,:;={}");
        char *svalue = strtok(NULL, " ,:;={}");
        // Reset allread variable
        allread = 0;
        // Sweep accross all parameters
        for (int i = 0; i < nprogpar; i++) {
            read_int_parameter(&progpar[i], key, svalue);
            // Add a value to allread if the parameter is already read
            if (progpar[i].read == 1) {
                allread = allread + 1;
            }
        }
        // If allread is equal to the number of parameters, break the loop
        if (allread == nprogpar) {
            break;
        }
    }
    // Check if all program input parameters were inserted
    check_for_int_missing_parameters(progpar, nprogpar);
}

void read_and_check_sys_parameters(FILE *input, char *buffer, int nsyspar, doublepar *syspar) {
    // Set the position to the beggining of the input file
    rewind(input);
    while(fgets(buffer, BUFSIZE, input)) {
        // Split line into the buffer with the following separators
        char *key = strtok(buffer, " ,:;={}");
        char *svalue = strtok(NULL, " ,:;={}");
        // Sweep accross all parameters
        for (int i = 0; i < nsyspar; i++) {
            read_double_parameter(&syspar[i], key, svalue);
        }
    }
    // Check if all program input parameters were inserted
    check_for_double_missing_parameters(syspar, nsyspar);
}

void read_and_check_nrms_param(FILE *input, char *buffer, intpar *nrmscalc, int limit) {
    // Set the position to the beggining of the input file
    rewind(input);
    while(fgets(buffer, BUFSIZE, input)) {
        // Split line into the buffer with the following separators
        char *key = strtok(buffer, " ,:;={}");
        char *svalue = strtok(NULL, " ,:;={}");
        read_int_parameter(nrmscalc, key, svalue);
    }    
    // Check if nrms is greater than dimension of the system
    check_for_invalid_int_input_value(nrmscalc, 1, limit, "dim (dimension / number of state variables of the system)", 0);
}

void read_and_check_nccalc_param(FILE *input, char *buffer, intpar *nccalc) {
    // Set the position to the beggining of the input file
    rewind(input);
    while(fgets(buffer, BUFSIZE, input)) {
        // Split line into the buffer with the following separators
        char *key = strtok(buffer, " ,:;={}");
        char *svalue = strtok(NULL, " ,:;={}");
        read_int_parameter(nccalc, key, svalue);
    }    
    // Check if nrms is greater than dimension of the system
    check_for_invalid_int_input_value(nccalc, 1, 0, "0", 1);
}

void read_and_check_ccalc_param(FILE *input, char *buffer, intpar *ccalc, int limit) {
    // Set the position to the beggining of the input file
    rewind(input);
    while(fgets(buffer, BUFSIZE, input)) {
        // Split line into the buffer with the following separators
        char *key = strtok(buffer, " ,:;={}");
        char *svalue = strtok(NULL, " ,:;={}");
        read_int_parameter(ccalc, key, svalue);
    }    
    // Check if nrms is greater than dimension of the system
    check_for_invalid_int_input_value(ccalc, 1, limit, "nccalc (number of custom calculations added to be performed)", 0);
}

void read_and_check_IC_params(char *ICstring, int dim, doublepar *ICpar) {
    /* Sweep ICstring to get values */
    // Declare svalue and get the first token "IC" to be discarted
    char *svalue = strtok(ICstring, " ,;:={}\n\r\t");
    // Get IC values
    int k = 0;
    while (svalue != NULL) {
        check_for_list_overflow(k, dim, "IC", "dim (dimension of the system)");
        // Get next value
        svalue = strtok(NULL, " ,;:={}\n\r\t");
        // Prevent Segmentation fault error if user forget to declare one or more parameters
        if (svalue != NULL) {
            read_double_parameter(&ICpar[k], ICpar[k].name, svalue);
        }    
        k = k + 1;
    }
    // Check if there is parameters missing in the input file
    check_for_double_missing_parameters(ICpar, dim);
}

void read_and_check_bifurc_params(char *bifstring, doublepar *bifpar) {
    /* Sweep bifstring to get values */
    // Declare svalue and get the first token "IC" to be discarted
    char *svalue = strtok(bifstring, " ,;:={}\n\r\t");
    // Get IC values
    int k = 0;
    while (svalue != NULL) {
        check_for_list_overflow(k, 3, "biflimits", "3");
        // Get next value
        svalue = strtok(NULL, " ,;:={}\n\r\t");
        // Prevent Segmentation fault error if user forget to declare one or more parameters
        if (svalue != NULL) {
            read_double_parameter(&bifpar[k], bifpar[k].name, svalue);
        }    
        k = k + 1;
    }
    // Check if there is parameters missing in the input file
    check_for_double_missing_parameters(bifpar, 3);
    check_for_invalid_double_input_value(&bifpar[2], 1, 0, "0", 1);
}

void read_and_check_rms_indexes(char *rmsstring, int nrms, intpar *rmspar, int limit) {
    /* Sweep rmsstring to get values */
    // Declare svalue and get the first token "IC" to be discarted
    char *svalue = strtok(rmsstring, " ,;:={}\n\r\t");
    // Get IC values
    int k = 0;
    while (svalue != NULL) {
        check_for_list_overflow(k, nrms, "rms", "nrms (number stated space variables submitted to rms calculations)");
        // Get next value
        svalue = strtok(NULL, " ,;:={}\n\r\t");
        // Prevent Segmentation fault error if user forget to declare one or more parameters
        if (svalue != NULL) {
            read_int_parameter(&rmspar[k], rmspar[k].name, svalue);
        }    
        k = k + 1;
    }
    // Check if there is parameters missing in the input file
    check_for_int_missing_parameters(rmspar, nrms);
    // Check if any element inserted in the rms list is bigger than the dimension of the system        
    check_for_invalid_int_input_value(rmspar, nrms, limit, "dim (dimension / number of state variables of the system)", 0);
}

void read_and_check_ccalc_indexes(char *ccalcstring, int nccalc, intpar *ccalc, char *listname, int limit, char *overflowdescription) {
    /* Sweep rmsstring to get values */
    // Declare svalue and get the first token "IC" to be discarted
    char *svalue = strtok(ccalcstring, " ,;:={}\n\r\t");
    // Get IC values
    int k = 0;
    while (svalue != NULL) {
        check_for_list_overflow(k, nccalc, listname, overflowdescription);
        // Get next value
        svalue = strtok(NULL, " ,;:={}\n\r\t");
        // Prevent Segmentation fault error if user forget to declare one or more parameters
        if (svalue != NULL) {
            read_int_parameter(&ccalc[k], ccalc[k].name, svalue);
        }    
        k = k + 1;
    }
    // Check if there is parameters missing in the input file
    check_for_int_missing_parameters(ccalc, nccalc);
    // Check if any element inserted in the rms list is bigger than the dimension of the system        
    check_for_invalid_int_input_value(ccalc, nccalc, limit, "nccalc (number of added customcalc to be performed).", 0);
}


int main () {
    // Open File and allocate memory for buffer
    char* name = "input.txt";
    FILE *input = fopen(name, "r");
    if (input == NULL) {
        // Return error if input does not exist 
        perror(name);
        exit(1);
    }
    char *buffer = malloc(BUFSIZE+1);
    /* ====================================================================== */
    // HANDLE PROGRAM PARAMETERS
    /* ====================================================================== */
    // Declare variable names
    char *progparnames[] = {"dim", "npar", "np", "ndiv", "trans", "maxper", 
                            "bifmode", "bifpar"};
    // Get the number of parameters in progparnames
    int nprogpar = sizeof(progparnames) / sizeof(progparnames[0]);
    // Declare Structure and Initialize it
    intpar progpar[nprogpar];
    initialize_intpar_struct(progpar, nprogpar, progparnames);
    // Read prog parameters and check if all were inserted in the input file
    read_and_check_prog_parameters(input, buffer, nprogpar, progpar);
    /* ====================================================================== */
    // HANDLE SYSTEM PARAMETERS
    /* ====================================================================== */
    // Search for the value that holds the number of parameters of the system
    int npar = search_for_vector_dimension(nprogpar, progpar, "npar");
    // Declare Struct to store IC parameters and initialize it 
    doublepar syspar[npar];
    initialize_doublepar_list_struct(syspar, npar, "par");
    // Read system parameters and check if all were inserted in the input file
    read_and_check_sys_parameters(input, buffer, npar, syspar);
    /* ====================================================================== */
    // HANDLE INITIAL CONDITIONS
    /* ====================================================================== */
    // Seartch for the value that holds the length of the IC vector
    int dim = search_for_vector_dimension(nprogpar, progpar, "dim");
    // Declare Struct and variables to store IC parameters and t0
    doublepar ICpar[dim];  // Initial conditions
    doublepar tpar[1];     // Initial Time
    char *tname[1] =  { "t0" };
    // Read file to store IC parameter list and check if it is inserted in the input file
    char* ICstring = read_and_check_list(input, buffer, "IC", "IC = {x[0], x[1], ...}");
    // Initialize IC structure
    initialize_doublepar_list_struct(ICpar, dim, "x");
    initialize_doublepar_struct(tpar, 1, tname);
    // Read IC and t0 parameter list and check if all parameters were inserted
    read_and_check_IC_params(ICstring, dim, ICpar);
    read_and_check_sys_parameters(input, buffer, 1, tpar);
    /* ====================================================================== */
    // HANDLE BIFURCATION PARAMETERS
    /* ====================================================================== */
    // Declare struct to store bifurcation limits and steps
    doublepar bifpar[3];       // Bifurcation limits and steps
    // Read file to store biflimits list and check if it is inserted in the input file
    char *bifstring = read_and_check_list(input, buffer, "biflimits", "biflimits = { bifpar[0], bifpar[1], bifpar[2] }");
    // Initialize bifpar structure
    initialize_doublepar_list_struct(bifpar, 3, "bifpar");
    // Read bif parameter list and check if all parameters in the list were inserted
    read_and_check_bifurc_params(bifstring, bifpar);
    /* ====================================================================== */
    // HANDLE RMS PARAMETERS
    /* ====================================================================== */
    // Declare struct to store rms list size
    intpar rmscalc;
    // Declare the name of the int variable that determine if there is rms calculations to be performed
    char *rmscalcname = "nrms";
    // Initialize struct of rmscalc
    initialize_intpar_struct(&rmscalc, 1, &rmscalcname);
    // Read and check 
    read_and_check_nrms_param(input, buffer, &rmscalc, dim);
    // Declare struct to store rms parameters
    intpar rmspar[rmscalc.value];
    // If nrms > 0, then read the rms parameters
    char *rmsstring;
    if (rmscalc.value > 0) {
        // Read file to store rms indexes list and check if it is inserted in the input file
        rmsstring = read_and_check_list(input, buffer, "rms", "rms = {rmsindex[0], rmsindex[1], ...}");
        // Initialize rms structure
        initialize_intpar_list_struct(rmspar, rmscalc.value, "rms");
        // Read rms index list and check if all parameters in the list were inserted
        read_and_check_rms_indexes(rmsstring, rmscalc.value, rmspar, dim);
    } 
    /* ====================================================================== */
    // HANDLE CUSTOMCALC PARAMETERS
    /* ====================================================================== */
    // Declare struct to store customcalc lists size
    intpar nccalc;
    intpar nendccalc;
    intpar nbodyccalc;
    // Declare the name of the int variable that determine if there is rms calculations to be performed
    char *nccalcname = "nccalc";
    char *endccalcname = "nccalc_";
    char *bodyccalcname = "nccalc*";
    // Initialize struct of customcalc lists
    initialize_intpar_struct(&nccalc, 1, &nccalcname);
    initialize_intpar_struct(&nendccalc, 1, &endccalcname);
    initialize_intpar_struct(&nbodyccalc, 1, &bodyccalcname);
    // Read and check 
    read_and_check_nccalc_param(input, buffer, &nccalc);
    if (nccalc.value > 0) {
        read_and_check_ccalc_param(input, buffer, &nendccalc, nccalc.value);
        read_and_check_ccalc_param(input, buffer, &nbodyccalc, nccalc.value);
    }
    // Declare struct to store customcalc and screencuscomcalc parameters
    intpar endccalcpar[nendccalc.value];
    intpar bodyccalcpar[nbodyccalc.value];
    // IF nccalc.value > 0, then perform the reading of customcalc parameters
    char *endccalcstring;
    char *bodyccalcstring;
    if (nccalc.value > 0) {
        // If endccalc.value > 0, then read the end type customcalc parameters
        if (nendccalc.value > 0) {
            // Read file to store customcalc indexes list and check if it is inserted in the input file
            endccalcstring = read_and_check_list(input, buffer, "ccalc_", "ccalc_ = {ccalc_[0], ccalc_[1], ...}");
            // Initialize endccalc structure
            initialize_intpar_list_struct(endccalcpar, nendccalc.value, "ccalc_");
            // Read ccalc index list and check if all parameters in the list were inserted
            read_and_check_ccalc_indexes(endccalcstring, nendccalc.value, endccalcpar, "calc_", nccalc.value, "nccalc_ (number of end type added customcalc to be printed).");
        }
        // If nbodyccalc.value > 0, then read the body type customcalc parameters
        if (nbodyccalc.value > 0) {
            // Read file to store customcalc indexes list and check if it is inserted in the input file
            bodyccalcstring = read_and_check_list(input, buffer, "ccalc*", "ccalc* = {ccalc*[0], ccalc*[1], ...}");
            // Initialize bodyccalc structure
            initialize_intpar_list_struct(bodyccalcpar, nbodyccalc.value, "ccalc*");
            // Read ccalc index list and check if all parameters in the list were inserted
            read_and_check_ccalc_indexes(bodyccalcstring, nbodyccalc.value, bodyccalcpar, "calc*", nccalc.value, "nccalc* (number of body type added customcalc to be printed).");
        }
    }
    
    // Print Values
    for (int j = 0; j < nprogpar; j++) {
        printf("progpar[%d] = %s = %d\n", j, progpar[j].name, progpar[j].value);
    }
    for (int j = 0; j < npar; j++) {
        printf("syspar[%d] = %s = %lf\n", j, syspar[j].name, syspar[j].value);
    }
    printf("tpar = %s = %lf\n", tpar[0].name, tpar[0].value);
    for (int i = 0; i < dim; i++) {
        printf("ICpar[%d] = %s = %lf\n", i,  ICpar[i].name, ICpar[i].value);
    }
    for (int i = 0; i < 3; i++) {
        printf("bifpar[%d] = %s = %lf\n", i,  bifpar[i].name, bifpar[i].value);
    }
    if (rmscalc.value > 0) { 
        printf("nrms = %s = %d\n", rmscalc.name, rmscalc.value);
        for (int i = 0; i < rmscalc.value; i ++) {
            printf("rmspar[%d] = %s = %d\n", i,  rmspar[i].name, rmspar[i].value);
        }
    }
    if (nccalc.value > 0) {
        printf("ncustomcalc = %s = %d\n", nccalc.name, nccalc.value);
        if (nendccalc.value > 0) {
            for (int i = 0; i < nendccalc.value; i++) {
                printf("endccalcpar = %s = %d\n", endccalcpar[i].name, endccalcpar[i].value);
            }
        }
        if (nbodyccalc.value > 0) {
            for (int i = 0; i < nbodyccalc.value; i++) {
                printf("bodyccalcpar = %s = %d\n", bodyccalcpar[i].name, bodyccalcpar[i].value);
            }
        }
    }
    // Free Memory
    free(ICstring); free(bifstring); 
    if (rmscalc.value > 0) {
        free(rmsstring);
    }
    if (nccalc.value > 0) {
        if (nendccalc.value > 0) {
            free(endccalcstring);
        }   
        if (nbodyccalc.value > 0) {
            free(bodyccalcstring);
        } 
    }
    free(buffer);
    return(0);
}
