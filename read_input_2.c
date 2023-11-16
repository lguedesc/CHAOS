#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "src/libs/msg.h"
#include "src/libs/iofiles.h"
#include "src/libs/basic.h"
#include "src/libs/defines.h"
#include "systems_placeholder.h"

/*
Steps:

1. Read contents of the file [done!]
2. Check all fields existance and for field duplicates [done!]
3. Check if the contents of the field are valid
   - If field if empty [done!]
   - If field is a string [done!]
   - If field is a number (int or double) [done!]
   - If field is a negative number [done!]
4. Handle system identification to read system parameters and IC
5. Handle lists 
6. Handle optional arguments
7. 'nrms' is not needed (read length of list and determine nrms (default length=0))

*/

/* ============================ STRUCTURES ============================= */

// Struct holding information about a single dynamical system
typedef struct {
    char *func_name;
    void (*func)(double, double);
    size_t dim;
    size_t npar;
    char *output_name;
    char *full_name;
    char *group;
    size_t angles;
    int *angle_indexes;
} ode_system;

// Function to find a system by name and return the function pointer
void (*find_ode_system(char *func_name, ode_system *system_list, size_t system_list_len))(double, double) {
    // Loop through system_list struct to find the corresponding name
    for (int i = 0; i < system_list_len; i++) {
        // Check if func_name is equal any of the system_list's func_names
        if (strcmp(func_name, system_list[i].func_name) == 0) {
            return system_list[i].func;
        }
    }
    // If it finds no corresponding name
    printf(COLOR_RED "System " COLOR_MAGENTA "'%s' " COLOR_RED "does not exist in the list of dynamical systems stored in this program version. Please add this new system or choose one of the following:\n", func_name);
    for (int i = 0; i < system_list_len; i++) {
        printf(COLOR_RED "   -> " COLOR_MAGENTA "'%s'\n", system_list[i].func_name);
    }
    printf("\n" RESET_STYLE);

    return NULL;
}

// Struct holding information about a single field (for each field)
typedef struct {
    char name[INPUT_BUFSIZE];
    void *value;
    bool duplicate;
    bool missing;
    bool valid;
    bool read;
} field;

field *init_field_struct(size_t length, char **field_names) {
    field *params = malloc(length * sizeof(*params));    
    // Check pointer args
    double_ptr_safety_check((void**)field_names, "char *field_names[] in 'static void init_field_struct(...)' function");    
    // Loop through all struct array spaces and initialize
    for (int i = 0; i < length; i++) {
        strcpy(params[i].name, field_names[i]);
        params[i].value = NULL;
        params[i].duplicate = false;
        params[i].missing = false;
        params[i].valid = false;
        params[i].read = false;
    }
    return params;
}

void free_field_struct(size_t length, field *params) {
    for (int i = 0; i < length; i++) {
        free(params[i].value);
    }
    free(params);
}

/* =========== FUNCTIONS TO HANDLE THE FILE AND ITS CONTENTS =========== */

FILE *open_and_check_file(char *filename) {
    // Open file
    FILE *file = fopen(filename, "r");
    // Check if the file opens successfully
    file_safety_check(file);
    return file;
}

size_t get_file_size(FILE *fp, char *filename) {
    // Start variable size
    size_t size = 0;
    // Check if file is NULL
    if (fp != NULL) {
        // Seek the end of the file
        fseek(fp, 0L, SEEK_END);
        // Return the current file position
        size = ftell(fp);
        // Go back to the beggining of the file
        rewind(fp);        
    }
    else {
        print_error("No such file or directory");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }

    return size;
}

char *read_input_file_and_store_contents(bool print_steps, bool print_file_contents) {
    // Get name/path of the file by input
    char *filename = get_input_filename();
    if (filename == NULL) {
        print_error("Error getting the filename!\n");
    }
    // Open and check file 
    FILE *file = open_and_check_file(filename);
    // Get file size
    size_t fsize = get_file_size(file, filename);
    printf(COLOR_YELLOW "=> " COLOR_GREEN "File " COLOR_MAGENTA "'%s' (%zu bytes) " COLOR_GREEN "opened successfully!\n" RESET_STYLE, filename, fsize);
    // Allocate memory for the contents of the file
    char *contents = malloc(fsize + 1);
    ptr_safety_check(contents, "'char *contents' in char *read_input_file_and_store_content() function");
    // Read the contents of the file into memory
    printf(COLOR_YELLOW "=> " RESET_STYLE "Reading file content...\n" RESET_STYLE);
    fread(contents, fsize + 1, 1, file);
    // Null-terminate the buffer to make it a valid C string
    contents[fsize] = '\0';
    printf(COLOR_YELLOW "=> " COLOR_GREEN "Content read successfully!\n" RESET_STYLE);
    // Print file contents
    if (print_file_contents == true) {
        print_warning("\nFILE CONTENTS:\n%s\n\n", contents);
    }
    // Close the file
    printf(COLOR_YELLOW "=> " RESET_STYLE "Closing file...\n" RESET_STYLE);
    fclose(file);
    printf(COLOR_YELLOW "=> " COLOR_GREEN "File closed successfully!\n" RESET_STYLE);


    return contents;
}

/* ================== FUNCTIONS TO GET INFORMATIONS ==================== */

char *get_line_with_specific_field(char* file_contents, char *field_name) {
    // Copy file_contents information to manipulate it without modifying the original string
    char *file_contents_copy = copy_pointer(file_contents, strlen(file_contents) + 1, sizeof(char));
    //print_warning("FILE_CONTENTS_COPY:\n%s\n", file_contents_copy);
    // Find the substring of file_contents_copy containing 'field_name' at its beggining
    char *substring = strstr(file_contents_copy, field_name);
    // Check if field_name really occurs in "substring"
    if (substring != NULL) {
        // Find the start of the line
        while (substring > file_contents_copy && *substring != '\n') {
            substring--;
        }
        // Move one character forward to exclude the newline character (if newline char exists)
        if (*substring == '\n') {
            substring++;
        }
        // Find the first occurrence of "\n" (end of the line)
        char *end_of_line = strstr(substring, "\n");
        if (end_of_line != NULL) {
            // Null-terminate the line to print only the part containing "field_name"
            *end_of_line = '\0';
        }
        // Make a copy of the line before freeing file_contents_copy
        char *line = copy_pointer(substring, strlen(substring) + 1, sizeof(char));
        // Free memory
        free(file_contents_copy);

        return line;
    } 
    else {
        // Free memory
        free(file_contents_copy);

        return NULL;
    }
}

/* ======================== SECURITY FUNCTIONS ========================= */
// Helper function to check if a character is a word boundary
static int is_word_boundary_char(char c) {
    // Add any characters that should be considered as word boundaries
    return isspace(c) || c == '=' || c == ':';
}

unsigned int count_field_occurrence(char *file_contents, char *field_name) {
    // Pointer to the first occurrence of the field name within file_contents
    const char *occurrence = file_contents;
    // Variable to keep track of occurrences
    int count = 0;
    // Length of the field name
    int field_len = strlen(field_name);
    // Iterate until no more occurrences are found
    while ((occurrence = strstr(occurrence, field_name)) != NULL) {
        // Check if the character before the occurrence is a word boundary
        bool is_word_boundary_before = (occurrence == file_contents) || is_word_boundary_char(*(occurrence - 1));
        // Check if the character after the occurrence is a word boundary
        bool is_word_boundary_after = is_word_boundary_char(*(occurrence + field_len)) || *(occurrence + field_len) == '\0';
        
        // Check if there are only blank spaces between the field_name and "=" sign
        const char *between_occurrence_and_equal = occurrence + field_len;
        while (isspace(*between_occurrence_and_equal)) {
            between_occurrence_and_equal++;
        }
        bool only_spaces_between = *between_occurrence_and_equal == '=';

        // If all conditions are true, it's a whole word match with only spaces between field_name and "="
        if (is_word_boundary_before && is_word_boundary_after && only_spaces_between) {
            // Move the pointer forward to avoid finding the same occurrence again
            occurrence += field_len;
            // Increment the count for each occurrence found
            count++;
        } else {
            // Move the pointer forward to search for the next occurrence
            occurrence++;
        }
    }

    return count;
}

bool is_field_missing(unsigned int n) {
    if (n == 0) {
        return true;
    }
    else {
        return false;
    }
}

bool is_field_duplicate(unsigned int n) {
    if (n > 1) {
        return true;
    }
    else {
        return false;
    }
}

void print_field_struct(field params, char *field_mode) {
    printf(COLOR_GREEN "name " COLOR_YELLOW "= " COLOR_MAGENTA "%16s" COLOR_WHITE " | " RESET_STYLE, params.name);
    printf(COLOR_GREEN "missing " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.missing PRINT_BOOL);
    printf(COLOR_GREEN "duplicate " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.duplicate PRINT_BOOL);
    printf(COLOR_GREEN "valid " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.valid PRINT_BOOL);
    printf(COLOR_GREEN "read " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.read PRINT_BOOL );
    printf(COLOR_GREEN "value " COLOR_YELLOW "= " COLOR_MAGENTA);
    if (strcmp(field_mode, "string") == 0) {
        printf("%s\n" RESET_STYLE, (char *)params.value);
    }
    else if (strcmp(field_mode, "int") == 0 || strcmp(field_mode, "int+") == 0) {
        printf("%d\n" RESET_STYLE, *(int *)params.value);
    }
    else if (strcmp(field_mode, "double") == 0 || strcmp(field_mode, "double+") == 0) {
        printf("%lf\n" RESET_STYLE, *(double *)params.value);
    }
    else if (strcmp(field_mode, "int_list") == 0 || strcmp(field_mode, "int_list+") == 0) {
        printf("--\n" RESET_STYLE);
    }
    else if (strcmp(field_mode, "double_list") == 0 || strcmp(field_mode, "double_list+") == 0) {
        printf("--\n" RESET_STYLE);
    }
    else {
        printf("--\n" RESET_STYLE);
    }

}

void print_field_mode(char *mode) {
    if (strcmp(mode, "string") == 0) {
        print_cyan("1. field_mode = %s\n", mode);
    }
    else if (strcmp(mode, "int") == 0) {
        print_cyan("2. field_mode = %s\n", mode);
    } 
    else if (strcmp(mode, "int+") == 0) {
        print_cyan("3. field_mode = %s\n", mode);
    }
    else if (strcmp(mode, "double") == 0) {
        print_cyan("4. field_mode = %s\n", mode);
    } 
    else if (strcmp(mode, "double+") == 0) {
        print_cyan("5. field_mode = %s\n", mode);
    } 
    else if (strcmp(mode, "int_list") == 0) {
        print_cyan("6. field_mode = %s\n", mode);
    }
    else if (strcmp(mode, "int_list+") == 0) {
        print_cyan("7. field_mode = %s\n", mode);
    }
    else if (strcmp(mode, "double_list") == 0) {
        print_cyan("8. field_mode = %s\n", mode);
    }
    else if (strcmp(mode, "double_list+") == 0) {
        print_cyan("9. field_mode = %s\n", mode);
    }
    else {
        print_cyan("field_mode = UNDEFINED\n");
    }
}

void print_check_messages(field param) {
    printf(COLOR_YELLOW "=> " COLOR_RED "The field " COLOR_MAGENTA "'%s' " COLOR_RED " not found in file contents.\n" RESET_STYLE, param.name);
}

bool is_field_value_missing(const char *field_line) {
    // Find occurrence position of "=" sign
    char *position = strchr(field_line, '=');
    // Check if equal sign is there
    if (position == NULL) {
        return true;
    }
    else {
        // neglect all white spaces after '='
        while (*(position + 1) == ' ') {
            position++;
        }
        // check if after all whitespaces it is a end of line character
        if (*(position + 1) == '\0') {
            return true;
        }
    }
    // If all conditions were not met
    return false;
}

char *get_field_value(char *field_line) {
    if (field_line != NULL) {
        // Copy the string to manipulate it without modifying the original
        char *field_line_copy = copy_pointer(field_line, strlen(field_line) + 1, sizeof(char));
        //print_blue("FIELD_LINE_COPY: %s\n", field_line_copy);
        // Tokenize the field_line_copy
        char *separators = "= ";
        char *field_value = strtok(field_line_copy, separators);
        if (field_value != NULL) {
            field_value = strtok(NULL, separators);
        }
        return field_value;
    }
    else {
        return NULL;
    }
    
}

bool is_field_value_contains_letter(const char *field_value) {
    // Check existance of field_value
    if (field_value != NULL) {
        // Lopp through field_value
        while (*field_value != '\0') {
            // Check if the character is a letter
            if (isalpha(*field_value)) { 
                return true;
            }
            // Move to the next character
            field_value++; 
        }
    }
    // If conditions weren't met, return false    
    return false;
}

bool is_string_number(char *str) {
    char *ptr;
    // Check if string is not NULL
    if (str != NULL) {
        // Try to convert string to double
        double number = strtod(str, &ptr);
        // Check if the pointers are equal (if no characters were converted, these pointers are equal)
        // and if the entire string was converted to ensure it is a whole number
        if (ptr == str || *ptr != '\n' && *ptr != '\0') {
            return false;
        }
        else {
            return true;
        }
    }
    else {
        return false;
    }
}

bool is_string_number_double(char *str_number) {
    // THIS FUNCTION SHOULD ONLY BE USED IF *str_number IS A NUMBER 
    // Check if string is not NULL
    if (str_number != NULL) {
        // Check if is a decimal point in the string
        if (strchr(str_number, '.') != NULL) {
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

bool is_string_number_positive(char *str_number) {
    // THIS FUNCTION SHOULD ONLY BE USED IF *str_number IS A NUMBER 
    // Check if str_number is not NULL
    if (str_number != NULL) {
        // Check if there is a negative sign in the first position of the string
        if (str_number[0] == '-') {
            return false;
        }
        else {
            return true;
        }
    }
    else {
        return false;
    }

}

bool is_field_value_valid(char *field_value, char *field_mode) {
    // Check if field_line is NULL or if field value is missing
    if (field_value == NULL) {
        return false;
    }
    // Check if field value is a number
    bool field_value_is_number = is_string_number(field_value);
    // Check mode to make compararisons
    if (strcmp(field_mode, "string") == 0) {
        // return true if field_value is not a number, and false if field_value is a number
        //print_blue("!field_value_is_number = %s\n", !field_value_is_number PRINT_BOOL);
        return !field_value_is_number;
    }
    else if (strcmp(field_mode, "int") == 0) {
        // returns true if field_value is a number AND field_value is not a double
        return field_value_is_number && !is_string_number_double(field_value);
    }
    else if (strcmp(field_mode, "int+") == 0) {
        // returns true if field_value is a number AND field_value is not a double AND field_value is positive
        return field_value_is_number && !is_string_number_double(field_value) && is_string_number_positive(field_value);
    }
    else if (strcmp(field_mode, "double") == 0) {
        // returns true if field_value is a number AND field_value is a double
        return field_value_is_number && is_string_number_double(field_value);
    }
    else if (strcmp(field_mode, "double+") == 0) {
        // returns true if field_value is a number AND field_value is a double AND field_value is positive
        return field_value_is_number && is_string_number_double(field_value) && is_string_number_positive(field_value);
    }
    else {
        print_debug("'char *field_mode' in 'is_field_value_valid()' function is invalid. Please check the code and enter one of the following: 'string', 'int', 'int+' 'double', 'double+', 'int_list', 'int_list+', 'double_list', 'double_list+'.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

void assign_field_value(field *param, char *field_value, char *field_mode) {
    // Assign field value based on mode:
    if (strcmp(field_mode, "string") == 0) {
        size_t len = strlen(field_value);
        param->value = malloc(len + 1);
        strcpy((*param).value, field_value);
    }
    else if (strcmp(field_mode, "int") == 0 || strcmp(field_mode, "int+") == 0) {
        int *tmp_value = malloc(sizeof *tmp_value);
        (*tmp_value) = atoi(field_value);
        param->value = tmp_value;
    }
    else if (strcmp(field_mode, "double") == 0 || strcmp(field_mode, "double+") == 0) {
        double *tmp_value = malloc(sizeof *tmp_value);
        (*tmp_value) = atof(field_value);
        param->value = tmp_value;
    }
    
    // Change read member to true if param.value is not null
    if ((*param).value != NULL) {
        (*param).read = true;
    }
}

/* ============= FUNCTIONS TO GET INFORMATION FROM SYSTEM ============== */

/* ==================== FULL CHECKS FOR PARAMETERS ===================== */

void check_fields(char *file_contents, size_t num_fields, field params[num_fields], char **field_names, ode_system *system_list, char *field_mode) {
    print_field_mode(field_mode);
    bool invalid = false;
    unsigned int n = 0;
    char *line = NULL;
    // Loop through all fields
    for (int i = 0; i < num_fields; i++) {
        // Copy field name to params[i].name member 
        if (strlen(field_names[i] + 1) <= INPUT_BUFSIZE) {
            strcpy(params[i].name, field_names[i]);
        }  
        else {
            print_debug("Field name '%s' with length greater than INPUT_BUFSIZE. Reduce the name.\n", field_names[i]);
        }
        // Count field occurrence
        n = count_field_occurrence(file_contents, field_names[i]);
        // Check for missing or duplicate field
        params[i].missing = is_field_missing(n);
        params[i].duplicate = is_field_duplicate(n);
        // Get the corresponding field line based on parameter (field) name and the field value
        line = get_line_with_specific_field(file_contents, params[i].name);
        char *field_value = get_field_value(line);
        //print_blue("|%s|\n", line);
        // Check if the field_value is valid based on field_mode
        params[i].valid = is_field_value_valid(field_value, field_mode);
        //print_blue("params[%d].valid = %s\n", i, params[i].valid PRINT_BOOL);
        // Check if all conditions are met to assign value to param struct
        if (params[i].missing == false && params[i].duplicate == false && params[i].valid == true) {
            assign_field_value(&params[i], field_value, field_mode);
        }

        print_field_struct(params[i], field_mode);
    }
    if (invalid == true) {
        print_error("Invalid declarations in the input file. Check previous messages and fix input file fields accordingly.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

/* ============================ MAIN FUNCTION ========================== */

int main(void) {
    /* ============================================================== */
    /* Define Systems and its additional information */
    /* ============================================================== */
    int lorenz_angle_indexes[] = {0, 1};
    ode_system systems[] = {
        // func_name, func(), dim, npar, output_name, full_name, group, angles, angle_indexes
        {"duffing", duffing, 2, 5, "duffing", "Duffing Oscillator", "HOS", 0, NULL},
        {"vanderpol", vanderpol, 2, 3, "vanderpol", "Van Der Pol Oscillator", "HOS", 0, NULL},
        {"lorenz", lorenz, 3, 4, "lorenz", "Lorenz System", "GNL", 2, lorenz_angle_indexes},
        {"bistable_EH", bistable_EH, 3, 8, "bistable_EH", "Polynomial Bistable Energy Harvester", "HOS", 0, NULL}                   
    };
    size_t systems_len = sizeof(systems)/sizeof(systems[0]);
    /* ============================================================== */
    /* Read File */
    /* ============================================================== */
    // Get the contents of the file
    char *fcontents = read_input_file_and_store_contents(true, false);
    /* ============================================================== */
    /* Declare Field Names */
    /* ============================================================== */
    char *system_field_names[] = {"group", "system"};
    size_t system_field_names_len = sizeof(system_field_names) / sizeof(system_field_names[0]);
    char *simulation_field_names[] = {"simulation", "integrator", "integration_mode"};
    size_t simulation_field_names_len = sizeof(simulation_field_names) / sizeof(simulation_field_names[0]);
    char *integrator_field_names[] = {"nP", "nDiv", "trans", "maxper"};
    size_t integrator_field_names_len = sizeof(integrator_field_names) / sizeof(integrator_field_names[0]);
    char *IC_field_names[] = {"IC", "t0"};
    size_t IC_field_names_len = sizeof(IC_field_names) / sizeof(IC_field_names[0]);
    char *RMScalc_field_names[] = {"nrms", "xRMS"};
    size_t RMScalc_fields_names_len = sizeof(RMScalc_field_names) / sizeof(RMScalc_field_names[0]);
    /* ============================================================== */
    /* Initialize Structs */
    /* ============================================================== */
    field *system_fields = init_field_struct(system_field_names_len, system_field_names);
    field *simulation_fields = init_field_struct(simulation_field_names_len, simulation_field_names);
    field *integrator_fields = init_field_struct(integrator_field_names_len, integrator_field_names);
    field *IC_fields = init_field_struct(IC_field_names_len, IC_field_names);
    field *RMScalc_fields = init_field_struct(RMScalc_fields_names_len, RMScalc_field_names);
    /* ============================================================== */
    /* Safety Checks in File Contents */
    /* ============================================================== */
    // Check for system_field_names
    printf(COLOR_YELLOW "=> " RESET_STYLE "Identifying system declared in file contents...\n" RESET_STYLE);
    check_fields(fcontents, system_field_names_len, system_fields, system_field_names, systems, "string");
    // Determine number of parameters and dimension of the system
    
    // Check for simulation_field_names
    printf(COLOR_YELLOW "=> " RESET_STYLE "Identifying simulation declared in file contents...\n" RESET_STYLE);
    check_fields(fcontents, simulation_field_names_len, simulation_fields, simulation_field_names, systems, "string");
    // Check for integrator_field_names
    printf(COLOR_YELLOW "=> " RESET_STYLE "Identifying remaining parameters declared in file contents...\n" RESET_STYLE);
    check_fields(fcontents, integrator_field_names_len, integrator_fields, integrator_field_names, systems, "int+");
    // Check for IC_fields_names
    //is_all_fields_present_once(fcontents, num_IC_fields, IC_fields, IC_field_names);
    // Check for RMScalc_field_names
    //is_all_fields_present_once(fcontents, num_RMScalc_fields, RMScalc_fields, RMScalc_field_names);
    /* ============================================================== */
    /* Test System determination */
    /* ============================================================== */
    void (*odesys)(double, double) = find_ode_system(system_fields[1].value, systems, systems_len);
    if (odesys != NULL) {
        num_integrator(odesys, 1.0, 2.0);
    }
    
    //strcpy(main_fields[0].name, "duffing");
    //printf("main_fields[0].name = %s\n", main_fields[0].name);
    /* ============================================================== */
    /* Free Memory */
    /* ============================================================== */
    free(fcontents);
    free_field_struct(system_field_names_len, system_fields);
    free_field_struct(simulation_field_names_len, simulation_fields); 
    free_field_struct(integrator_field_names_len, integrator_fields); 
    free_field_struct(IC_field_names_len, IC_fields); 
    free_field_struct(RMScalc_fields_names_len, RMScalc_fields); 
}