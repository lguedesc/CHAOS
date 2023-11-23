#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "src/libs/msg.h"
#include "src/libs/iofiles.h"
#include "src/libs/basic.h"
#include "src/libs/defines.h"
#include "systems_placeholder.h"

#define MAX_INPUT_FIELD_VALUE_LEN 200
#define MAX_INPUT_FIELD_LEN 20
#define MAX_NUMBER_LEN_INSIDE_BRACKET 5

/* ============================ STRUCTURES ============================= */

// Struct holding information about a single dynamical system
typedef struct {
    char *func_name;
    size_t dim;
    size_t npar;
    char *full_name;
    char *group;
    void (*func)(int, double *, double, double *, double *);
    void (*customcalc_func)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int);
    size_t angles;
    int *angle_indexes;
} ode_system;

// Function to find a system by name and return the function pointer
void (*find_ode_system(char *func_name, ode_system *system_list, size_t system_list_len))(int, double *, double, double *, double *) {
    // Loop through system_list struct to find the corresponding name
    for (int i = 0; i < system_list_len; i++) {
        // Check if func_name is equal any of the system_list's func_names
        if (strcmp(func_name, system_list[i].func_name) == 0) {
            return system_list[i].func;
        }
    }
    // If it finds no corresponding name
    printf(COLOR_YELLOW "=> " COLOR_RED "System " COLOR_MAGENTA "'%s' " COLOR_RED "does not exist in the list of dynamical systems stored in this program version. Please add this new system or choose one of the following:\n", func_name);
    for (int i = 0; i < system_list_len; i++) {
        printf(COLOR_RED "   -> " COLOR_MAGENTA "'%s'\n", system_list[i].func_name);
    }
    printf("\n" RESET_STYLE);

    print_exit_prog();
    exit(EXIT_FAILURE);

    //return NULL;
}

void (*find_customcalc_func(char *func_name, ode_system *system_list, size_t system_list_len))(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int) {
    // Loop through system_list struct to find the corresponding name
    for (int i = 0; i < system_list_len; i++) {
        // Check if func_name is equal any of the system_list's func_names
        if (strcmp(func_name, system_list[i].func_name) == 0) {
            return system_list[i].customcalc_func;
        }
    }
    // If it finds no corresponding name
    print_debug(COLOR_RED "Customcalc function of system " COLOR_MAGENTA "'%s' " COLOR_RED "does not exist in the list of dynamical systems stored in this program version. Please check definition in source code.\n", func_name);

    print_exit_prog();
    exit(EXIT_FAILURE);
}

// Struct holding information about a single field (for each field)
typedef struct {
    char name[INPUT_BUFSIZE];
    void *value;
    bool duplicate;
    bool missing;
    bool optional;
    bool valid;
    bool errors;
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
        params[i].optional = false;
        params[i].valid = false;
        params[i].errors = false;
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

// Define a struct to hold a pointer to the fields and their length
typedef struct {
    field *fields;
    size_t len;
} field_group;

/* ========================= PRINT FUNCTIONS  ========================== */

void print_field_struct(field params, char *field_mode) {
    printf(COLOR_GREEN "name " COLOR_YELLOW "= " COLOR_MAGENTA "%16s" COLOR_WHITE " | " RESET_STYLE, params.name);
    printf(COLOR_GREEN "missing " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.missing PRINT_BOOL);
    printf(COLOR_GREEN "duplicate " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.duplicate PRINT_BOOL);
    printf(COLOR_GREEN "optional " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.optional PRINT_BOOL);
    printf(COLOR_GREEN "valid " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.valid PRINT_BOOL);
    printf(COLOR_GREEN "errors " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.errors PRINT_BOOL);
    printf(COLOR_GREEN "read " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.read PRINT_BOOL );
    
    printf(COLOR_GREEN "value " COLOR_YELLOW "= " COLOR_MAGENTA);
    if (params.value != NULL) {
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
    } else {
        printf("NULL\n" RESET_STYLE);
    }
}

void print_field_list_struct(field params, size_t n_cells, char *field_mode) {
    printf(COLOR_GREEN "name " COLOR_YELLOW "= " COLOR_MAGENTA "%16s" COLOR_WHITE " | " RESET_STYLE, params.name);
    printf(COLOR_GREEN "missing " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.missing PRINT_BOOL);
    printf(COLOR_GREEN "duplicate " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.duplicate PRINT_BOOL);
    printf(COLOR_GREEN "optional " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.optional PRINT_BOOL);
    printf(COLOR_GREEN "valid " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.valid PRINT_BOOL);
    printf(COLOR_GREEN "errors " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.errors PRINT_BOOL);
    printf(COLOR_GREEN "read " COLOR_YELLOW "= " COLOR_MAGENTA "%5s" COLOR_WHITE " | " RESET_STYLE, params.read PRINT_BOOL );
    printf(COLOR_GREEN "value " COLOR_YELLOW "= " COLOR_MAGENTA);
    if (params.value != NULL) {
        if (strcmp(field_mode, "int_list") == 0 || strcmp(field_mode, "int_list+") == 0) {
            printf("{ ");
            for (int i = 0; i < n_cells; i++) {
                if (i == n_cells - 1) {
                    printf("%d", ((int *)params.value)[i]);    
                }
                else {
                    printf("%d, ", ((int *)params.value)[i]);
                }
            }
            printf(" }\n" RESET_STYLE);
        }
        else if (strcmp(field_mode, "double_list") == 0 || strcmp(field_mode, "double_list+") == 0) {
            printf("{ ");
            for (int i = 0; i < n_cells; i++) {
                if (i == n_cells - 1) {
                    printf("%g", ((double *)params.value)[i]);    
                }
                else {
                    printf("%g, ", ((double *)params.value)[i]);
                }
            }
            printf(" }\n" RESET_STYLE);
        }
        else {
            printf("--\n" RESET_STYLE);
        }
    } else {
        printf("NULL\n" RESET_STYLE);
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

char *get_field_content(char* file_contents, char *field_name, char *field_mode) {
    // Security
    if (file_contents == NULL || field_name == NULL || field_mode == NULL) {
        return NULL;
    }
    // Copy file_contents information to manipulate it without modifying the original string
    char *file_contents_copy = copy_pointer(file_contents, strlen(file_contents) + 1, sizeof(char));
    //print_warning("FILE_CONTENTS_COPY:\n%s\n", file_contents_copy);
    // Find the substring of file_contents_copy containing 'field_name' at its beggining
    char *substring = strstr(file_contents_copy, field_name);
    // Check if field_name really occurs in "substring"
    if (substring != NULL) {
        // Find the start of the field content
        while (substring > file_contents_copy && *substring != '\n') {
            substring--;
        }
        // Move one character forward to exclude the newline character (if newline char exists)
        if (*substring == '\n') {
            substring++;
        }
        // Check field mode to determine end of content character
        if (strcmp(field_mode, "int_list") == 0 || strcmp(field_mode, "int_list+") == 0 || strcmp(field_mode, "double_list") == 0 || strcmp(field_mode, "double_list+") == 0){
            // Find the first occurrence of "}" (end of the content)
            char *end_of_content = strchr(substring, '}');
            if (end_of_content != NULL) {
                // Check if there is a second "=" sign between the beginning of the substring and the first occurrence of the character '}'
                char *first_equal_sign = strchr(substring, '=');
                char *second_equal_sign = strchr(first_equal_sign + 1, '=');
                if (second_equal_sign != NULL && second_equal_sign < end_of_content) {
                    // Move one character forward
                    first_equal_sign++;
                    // Null-terminate the content to print only the part containing the field name + =
                    *first_equal_sign = '\0';
                } else {
                    // Move one character forward
                    end_of_content++;
                    // Null-terminate the content to print only to part containing the field content
                    *end_of_content = '\0';
                }
            }
        }
        else {
            // Find the first occurrence of "\n" (end of the content)
            char *end_of_content = strstr(substring, "\n");
            if (end_of_content != NULL) {
                // Null-terminate the line to print only the part containing "field_name"
                *end_of_content = '\0';
            }
        }
        // Make a copy of the substring before freeing file_contents_copy
        char *field_content = copy_pointer(substring, strlen(substring) + 1, sizeof(char));
        // Free memory
        free(file_contents_copy);

        return field_content;
    } 
    else {
        // Free memory
        free(file_contents_copy);

        return NULL;
    }
}

char *get_field_value(char *field_contents) {
    if (field_contents != NULL) {
        // Copy the string to manipulate it without modifying the original
        char *field_contents_copy = copy_pointer(field_contents, strlen(field_contents) + 1, sizeof(char));
        //print_blue("FIELD_CONTENTS_COPY: %s\n", field_contents_copy);
        // Tokenize the field_contents_copy
        char *separators = "= ";
        char *field_value = strtok(field_contents_copy, separators);
        if (field_value != NULL) {
            field_value = strtok(NULL, separators);
        }
        return field_value;
    }
    else {
        return NULL;
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
    else {
        print_debug("Invalid 'char *field_mode' in ' void assign_field_value(...)' function.\n");
    }
    // Change read member to true if param.value is not null
    if ((*param).value != NULL) {
        (*param).read = true;
    }
}

void assign_field_list_values(field *param, char **field_values_array, size_t n_cells, char *field_mode) {
    // Assign field value based on mode
    if (strcmp(field_mode, "int_list") == 0 || strcmp(field_mode, "int_list+") == 0) {
        int *tmp_values = malloc(n_cells * sizeof *tmp_values);
        // Iterate over the values 
        for (int i = 0; i < n_cells; i++) {
            tmp_values[i] = atoi(field_values_array[i]);
        }
        // Assign values to struct
        param->value = tmp_values;
    }
    else if (strcmp(field_mode, "double_list") == 0 || strcmp(field_mode, "double_list+") == 0) {
        double *tmp_values = malloc(n_cells * sizeof *tmp_values);
        // Iterate over the values 
        for (int i = 0; i < n_cells; i++) {
            tmp_values[i] = atof(field_values_array[i]);
        }
        // Assign values to struct
        param->value = tmp_values;
    }
    else {
        print_debug("Invalid 'char *field_mode' in ' void assign_field_list_values(...)' function.\n");
    }
    // Change read member to true if param.value is not null
    if (param->value != NULL) {
        (*param).read = true;
    }
}

char *get_field_list_value(char *field_contents) {
    // Check if field contents is NULL
    if (field_contents == NULL) {
        return NULL;
    }
    // Copy the string to manipulate it without modifying the original
    char *field_contents_copy = copy_pointer(field_contents, strlen(field_contents) + 1, sizeof(char));
    //print_blue("FIELD_CONTENTS_COPY: %s\n", field_contents_copy);
    // Get substring containing "{" at its beggining
    char *substring = strchr(field_contents_copy, '{');
    // Check if substring is not NULL
    if (substring != NULL) {
        // Find the first occurrence of } to find the end of the field value
        char *end_of_content = strchr(substring, '}');
        // Check if end of content is not NULL
        if (end_of_content != NULL) {
            // Move one character forward
            end_of_content++;
            // Null-terminate the content to print only to part containing the field content
            *end_of_content = '\0';
        }
        // Make a copy of the substring before freeing field_contents_copy
        char *field_list_value = copy_pointer(substring, strlen(substring) + 1, sizeof(char));
        // Free memory
        free(field_contents_copy);
        return field_list_value;
    }
    else {
        // Free memory
        free(field_contents_copy);
        return NULL;
    }
}

size_t get_number_of_list_cells_declared(char *field_value) {
    // Declare counter
    size_t count = 0;
    if (field_value == NULL) {
        return count;
    }
    else {
        // Iterate through each character in the field_value
        while (field_value != NULL && *field_value != '\0') {
            // Check if the current character matches the separator ","
            if (*field_value == ',') {
                count++;
            }
            // Move the to the next character in field_value
            field_value++;
        }
        // The number of cells is count (number of ",") + 1
        return count + 1;
    }
}

char **get_each_element_of_list(char *list, size_t n_cells, char *list_name) {
    // Check if there is any number of elements to allocate memory
    if (n_cells == 0 || list == NULL || list_name == NULL) {
        return NULL;
    } 
    else {
        // Make a copy of the list to modify it without changing the original
        char *list_copy = copy_pointer(list, strlen(list) + 1, sizeof(char));
        // Create an array of strings to store each element of the list
        char **list_elements = alloc_string_array(n_cells, MAX_INPUT_FIELD_VALUE_LEN);
        // Get n_cells values based on a separator
        char *token = NULL;
        char *separators = "{,} ";
        for (int i = 0; i < n_cells; i++) {
            if (i == 0) {
                token = strtok(list_copy, separators);
            } else {
                token = strtok(NULL, separators);
            }
            if (token == NULL) {
                strcpy(list_elements[i], " ");
            } else {
                if (strlen(token) > MAX_INPUT_FIELD_VALUE_LEN - 1) {
                    print_error("Length of the elements in '%s' list exceeds maximum length allowed. Each cell of the list must have less than %d characters.\n", list_name, MAX_INPUT_FIELD_VALUE_LEN);
                    print_exit_prog();
                    exit(EXIT_FAILURE);
                } else {
                    strcpy(list_elements[i], token);
                }
            }
        }
        return list_elements;
    }
}

char **get_names_of_table(size_t npar, char *name_of_cell) {
    // Check if there is any number of elements to allocate memory
    if (npar == 0 || name_of_cell == NULL) {
        return NULL;
    }
    else {
        char **table_element_names = NULL;
        if (strlen(name_of_cell) > MAX_INPUT_FIELD_VALUE_LEN - MAX_NUMBER_LEN_INSIDE_BRACKET - 2) {
            print_debug("String in 'char *name_of_cell' too long in char **get_names_of_table(...) function\n");   
            print_exit_prog();
            exit(EXIT_FAILURE);
        } 
        else {
            // Create an array of strings to store each element of the table
            table_element_names = alloc_string_array(npar, MAX_INPUT_FIELD_VALUE_LEN);
            // Assign string names
            for (int i = 0; i < npar; i++) {
                sprintf(table_element_names[i], "%s[%d]", name_of_cell, i);
            }
        }
        return table_element_names;
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

bool is_list_overflow(size_t n_cells, size_t expected_n_cells) {
    // Check if n_elements is greater than n_cells
    if (n_cells > expected_n_cells) {
        return true;
    }
    else {
        return false;
    }
}

bool is_list_underflow (size_t n_cells, size_t expected_n_cells) {
    // Check if expected_n_elements is lesser than n_cells
    if (n_cells < expected_n_cells) {
        return true;
    } 
    else {
        return false;
    }
}

bool is_each_element_of_list_a_number(char **array_with_all_elements, size_t n_cells) {
    // Iterate over element list to check if the elements are numbers or not
    for (int i = 0; i < n_cells; i++) {
        if (is_string_number(array_with_all_elements[i]) == false) {
            // If it is not number, return false
            return false;
        }
    }
    // Else, return true
    return true;
}

bool is_each_element_of_list_a_double_number(char **array_with_all_elements, size_t n_cells) {
    /* WARNING: THIS FUNCTION MUST BE USED WITH AN ARRAY OF STRINGS THAT ARE NUMBERS */
    // Iterate over element list to check if the elements are double numbers or not
    for (int i = 0; i < n_cells; i++) {
        if (is_string_number_double(array_with_all_elements[i]) == false) {
            // If it is not double, return false
            return false;
        }
    }
    // Else, return true
    return true;
}

bool is_each_element_of_list_an_int_number(char **array_with_all_elements, size_t n_cells) {
    /* WARNING: THIS FUNCTION MUST BE USED WITH AN ARRAY OF STRINGS THAT ARE NUMBERS */
    // Iterate over element list to check if the elements are int numbers or not
    for (int i = 0; i < n_cells; i++) {
        if (is_string_number_double(array_with_all_elements[i]) == true) {
            // If it is double, return false
            return false;
        }
    }
    // Else, return true
    return true;
}

bool is_each_element_of_list_a_positive_number(char **array_with_all_elements, size_t n_cells) {
    /* WARNING: THIS FUNCTION MUST BE USED WITH AN ARRAY OF STRINGS THAT ARE NUMBERS */
    // Iterate over element list to check if the elements are all positive or not
    for (int i = 0; i < n_cells; i++) {
        if (is_string_number_positive(array_with_all_elements[i]) == false) {
            // If it is not number, return false
            return false;
        }
    }
    // Else, return true
    return true;
}

bool is_limit_crossed(char *str_number, size_t *limit, char *field_mode) {
    /* ONLY USE THIS FUNCTION IF str_number IS INDEED A NUMBER */
    // Check if str_number is null
    if (str_number != NULL && limit != NULL) {
        // Check variable type
        if (strcmp(field_mode, "int") == 0 || strcmp(field_mode, "int+") == 0) {
            int number = atoi(str_number);
            if (number > (int)(*limit)) {
                return true;
            }
        }
        else if (strcmp(field_mode, "double") == 0 || strcmp(field_mode, "double+") == 0) {
            double number = atof(str_number);
            if (number > (double)(*limit)) {
                return true;
            }
        }
        else {
            print_debug("Invalid 'char *var_type' in 'bool is_limit_crossed(...)' function.\n");
            print_exit_prog();
            exit(EXIT_FAILURE);
        }
    } 
    
    return false;
}

bool does_list_contain_duplicates(char **field_values_array, size_t n_cells) {
    // Iterate over the field_values string array and compare all of them
    for (int i = 0; i < n_cells; i++) {
        for (int j = 0; j < n_cells; j++) {
            if(strcmp(field_values_array[i], field_values_array[j]) == 0 && i != j) {
                return true;
            }
        }
    }
    return false;
}

bool does_any_element_of_list_cross_a_limit(char **field_values_array, size_t n_cells, size_t *limit, char *field_mode) {
    // Check pointers before continuing
    if (field_values_array != NULL && *field_values_array != NULL && limit != NULL && field_mode != NULL) {
        // Check field_mode before continuing
        if (strcmp(field_mode, "int_list") == 0 || strcmp(field_mode, "int_list+") == 0) {
            int number;  
            // Iterate over the field_values string array to check if any value crossed the limit
            for (int i = 0; i < n_cells; i++) { 
                number = atoi(field_values_array[i]);
                if (number > (int)(*limit)) {
                    return true;
                }    
            }
        }
        else if (strcmp(field_mode, "double_list") == 0 || strcmp(field_mode, "double_list+") == 0) {
            double number;  
            // Iterate over the field_values string array to check if any value crossed the limit
            for (int i = 0; i < n_cells; i++) { 
                number = atof(field_values_array[i]);
                if (number > (double)(*limit)) {
                    return true;
                }    
            }
        }
        else {
            print_debug("Invalid 'char *field_mode' in 'bool does_any_element_of_list_cross_a_limit(...)' function.\n");
            print_exit_prog();
            exit(EXIT_FAILURE);
        }
    } 

    return false;
}

char *add_error_message(char *error_messages, const char *formatted_message, ...) {
    va_list args;
    va_start(args, formatted_message);

    // Determine the length of the formatted string
    va_list tmp_args;
    va_copy(tmp_args, args);
    int len = vsnprintf(NULL, 0, formatted_message, tmp_args);
    va_end(tmp_args);

    // Realloc error messages array based on its previous size and the length of the formatted string
    if (error_messages == NULL) {
        error_messages = malloc((len + 1) * sizeof(char));
    } 
    else {
        error_messages = realloc(error_messages, (strlen(error_messages) + len + 1) * sizeof(char));
    }

    // Append the formatted string to error_messages
    vsprintf(error_messages + strlen(error_messages), formatted_message, args);

    va_end(args);

    return error_messages;
}

void error_messages_for_invalid_field(bool is_number, bool is_double, bool limit_crossed, bool is_positive, size_t *limit, char **error_messages, char *field_mode) {
    // Add messages if needed
    if (strcmp(field_mode, "string") == 0) {
        // Add error messages if needed
        if (is_number == true) {
            (*error_messages) = add_error_message((*error_messages), "   - Field value is a number, it must be a string.\n");
        }
    }
    else if (strcmp(field_mode, "int") == 0 || strcmp(field_mode, "int+") == 0 ) {
        // Add error messages if needed
        if (is_number == false) {
            (*error_messages) = add_error_message((*error_messages), "   - Field value must be a number.\n");
        }
        if (is_double == true) {
            (*error_messages) = add_error_message((*error_messages), "   - Field value must be an integer number.\n");
        }
        if (limit_crossed == true) {
            (*error_messages) = add_error_message((*error_messages), "   - Field can't have values above %zu.\n", (*limit));
        }
        if (strcmp(field_mode, "int+") == 0 && is_positive == false) {
            (*error_messages) = add_error_message((*error_messages), "   - Field must be positive.\n");
        }
    }
    else if (strcmp(field_mode, "double") == 0 || strcmp(field_mode, "double+") == 0 ) {
        // Add error messages if needed
        if (is_number == false) {
            (*error_messages) = add_error_message((*error_messages), "   - Field value must be a number.\n");
        }
        if (is_double == false) {
            (*error_messages) = add_error_message((*error_messages), "   - Field value must be an floating point number (i.e.: 5.0 or 4.5).\n");
        }
        if (limit_crossed == true) {
            (*error_messages) = add_error_message((*error_messages), "   - Field can't have values above %zu.\n", (*limit));
        }
        if (strcmp(field_mode, "double+") == 0 && is_positive == false) {
            (*error_messages) = add_error_message((*error_messages), "   - Field must be positive.\n");
        }
    }
}

void error_messages_for_invalid_list(bool is_number_array, bool overflow, bool underflow, bool is_all_double, bool is_all_int, bool allow_cell_duplicates, bool contains_duplicates, bool is_each_element_positive, size_t expected_n_cells, bool any_element_cross_limit, size_t *cell_limit_value, char **error_messages, char *field_mode) {
    // Add error messages if needed 
    if (strcmp(field_mode, "int_list") == 0 || strcmp(field_mode, "int_list+") == 0 ) {
        if (is_number_array == false) {
            (*error_messages) = add_error_message((*error_messages), "   - List must have numbers in all elements %zu.\n");
        }
        if (overflow == true || underflow == true) {
            (*error_messages) = add_error_message((*error_messages), "   - List must have %zu elements.\n", expected_n_cells);
        }
        if (is_all_int == false) {
            (*error_messages) = add_error_message((*error_messages), "   - All elements of list must be integer numbers (i.e.: 1 or 5).\n");
        }
        if (cell_limit_value != NULL && any_element_cross_limit == true) {
            (*error_messages) = add_error_message((*error_messages), "   - Each element of list must have a maximum value of %zu.\n", (*cell_limit_value));
        }
        if (strcmp(field_mode, "int_list+") == 0 && is_each_element_positive == false) {
            (*error_messages) = add_error_message((*error_messages), "   - All elements of list must be positive.\n");
        }
        if (allow_cell_duplicates == false && contains_duplicates == true) {
            (*error_messages) = add_error_message((*error_messages), "   - List contains duplicates.\n");
        }
    }
    else if (strcmp(field_mode, "double_list") == 0 || strcmp(field_mode, "double_list+") == 0 ) {
        if (is_number_array == false) {
            (*error_messages) = add_error_message((*error_messages), "   - List must have numbers in all elements.\n");
        }
        if (overflow == true || underflow == true) {
            (*error_messages) = add_error_message((*error_messages), "   - List must have %zu elements.\n", expected_n_cells);
        }
        if (is_all_double == false) {
            (*error_messages) = add_error_message((*error_messages), "   - All elements of list must be double numbers (i.e.: 5.0 or 4.5).\n");
        }
        if (strcmp(field_mode, "double_list+") == 0 && is_each_element_positive == false) {
            (*error_messages) = add_error_message((*error_messages), "   - All elements of list must be positive.\n");
        }
        if (allow_cell_duplicates == false && contains_duplicates == true) {
            (*error_messages) = add_error_message((*error_messages), "   - List contains duplicates.\n");
        }
    }

}

bool is_field_value_valid(char *field_value, size_t *limit, char *field_mode, char **error_messages) {
    // Check if field_value is NULL or if field value is missing
    if (field_value == NULL) {
        print_debug("Failed to read 'char **field_value' in 'bool is_field_valid(...)' function.\n");
        return false;
    }
    // Check if field value is a number and a double
    bool is_number = is_string_number(field_value);
    bool is_double = is_string_number_double(field_value);
    bool limit_crossed = is_limit_crossed(field_value, limit, field_mode);
    bool is_positive = false;
    // Check mode to make compararisons
    if (strcmp(field_mode, "string") == 0) {
        // Add error messages if needed
        error_messages_for_invalid_field(is_number, is_double, limit_crossed, is_positive, &(*limit), error_messages, field_mode);
        // return true if field_value is not a number, and false if field_value is a number
        return !is_number;
    }
    else if (strcmp(field_mode, "int") == 0) {
        // returns true if field_value is a number AND field_value is not a double AND field_valie not crossed the limit
        return is_number && !is_double && !limit_crossed;
    }
    else if (strcmp(field_mode, "int+") == 0) {
        // Check if it is positive
        is_positive = is_string_number_positive(field_value);
        // Add error messages if needed
        error_messages_for_invalid_field(is_number, is_double, limit_crossed, is_positive, &(*limit), error_messages, field_mode);
        // returns true if field_value is a number AND field_value is not a double AND field_value is positive AND field_valie not crossed the limit
        return is_number && !is_double && is_positive && !limit_crossed;;
    }
    else if (strcmp(field_mode, "double") == 0) {
        // Add error messages if needed
        error_messages_for_invalid_field(is_number, is_double, limit_crossed, is_positive, &(*limit), error_messages, field_mode);
        // returns true if field_value is a number AND field_value is a double AND field_valie not crossed the limit
        return is_number && is_double && !limit_crossed;;
    }
    else if (strcmp(field_mode, "double+") == 0) {
        // Check if it is positive
        is_positive = is_string_number_positive(field_value);
        // Add error messages if needed
        error_messages_for_invalid_field(is_number, is_double, limit_crossed, is_positive, &(*limit), error_messages, field_mode);
        // returns true if field_value is a number AND field_value is a double AND field_value is positive AND field_valie not crossed the limit
        return is_number && is_double && is_positive && !limit_crossed;;
    }
    else {
        print_debug("'char *field_mode' in 'is_field_value_valid()' function is invalid. Please check the code and enter one of the following: 'string', 'int', 'int+' 'double', 'double+'.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

bool is_field_list_valid(char **field_values_array, size_t n_cells, size_t expected_n_cells, bool allow_cell_duplicates, size_t *cell_limit_value, char **error_messages, char *field_mode) {
    // Check if field_values_array is NULL
    if (field_values_array == NULL || *field_values_array == NULL) {
        print_debug("Failed to read 'char **field_values_array' in 'bool is_field_valid_list(...)' function.\n");
        return false;
    }
    // Check if each element of the field_values_array is a number
    bool is_number_array = is_each_element_of_list_a_number(field_values_array, n_cells);
    bool overflow = is_list_overflow(n_cells, expected_n_cells);
    bool underflow = is_list_underflow(n_cells, expected_n_cells);
    bool is_all_double = is_each_element_of_list_a_double_number(field_values_array, n_cells);
    bool any_element_cross_limit = does_any_element_of_list_cross_a_limit(field_values_array, n_cells, cell_limit_value, field_mode);
    bool is_each_element_positive = false;
    bool contains_duplicates = does_list_contain_duplicates(field_values_array, n_cells);
    bool is_all_int = is_each_element_of_list_an_int_number(field_values_array, n_cells);
    // Check for duplicates
    if (allow_cell_duplicates == false && contains_duplicates == true) {
        // Add error message if needed
        error_messages_for_invalid_list(is_number_array, overflow, underflow, is_all_double, is_all_int, allow_cell_duplicates, contains_duplicates, is_each_element_positive, expected_n_cells, any_element_cross_limit, cell_limit_value, error_messages, field_mode);
        return false;
    }
    // Check mode and make comparisons
    if (strcmp(field_mode, "int_list") == 0) {
        // Add error message if needed
        error_messages_for_invalid_list(is_number_array, overflow, underflow, is_all_double, is_all_int, allow_cell_duplicates, contains_duplicates, is_each_element_positive, expected_n_cells, any_element_cross_limit, cell_limit_value, error_messages, field_mode);
        // Return if all elements are numbers, if the list is not overflow or underflow, if all elements of the list are not doubles, and if limit is not crossed
        return is_number_array && !overflow && !underflow && is_all_int && !any_element_cross_limit;
    }
    else if (strcmp(field_mode, "int_list+") == 0) {
        is_each_element_positive = is_each_element_of_list_a_positive_number(field_values_array, n_cells);
        // Add error message if needed
        error_messages_for_invalid_list(is_number_array, overflow, underflow, is_all_double, is_all_int, allow_cell_duplicates, contains_duplicates, is_each_element_positive, expected_n_cells, any_element_cross_limit, cell_limit_value, error_messages, field_mode);
        // Return if all elements are numbers, if the list is not overflow or underflow, and if all elements of the list are not doubles, and if all numbers are positive, and if limit is not crossed
        return is_number_array && !overflow && !underflow && is_all_int && is_each_element_positive && !any_element_cross_limit;
    }
    else if (strcmp(field_mode, "double_list") == 0) {
        // Add error message if needed
        error_messages_for_invalid_list(is_number_array, overflow, underflow, is_all_double, is_all_int, allow_cell_duplicates, contains_duplicates, is_each_element_positive, expected_n_cells, any_element_cross_limit, cell_limit_value, error_messages, field_mode);
        // Return if all elements are numbers, if the list is not overflow or underflow, and if all elements of the list are doubles, , and if limit is not crossed
        return is_number_array && !overflow && !underflow && is_all_double && !any_element_cross_limit;;
    }
    else if (strcmp(field_mode, "double_list+") == 0) {
        is_each_element_positive = is_each_element_of_list_a_positive_number(field_values_array, n_cells);
        // Add error message if needed
        error_messages_for_invalid_list(is_number_array, overflow, underflow, is_all_double, is_all_int, allow_cell_duplicates, contains_duplicates, is_each_element_positive, expected_n_cells, any_element_cross_limit, cell_limit_value, error_messages, field_mode);
        // Return if all elements are numbers, if the list is not overflow or underflow, and if all elements of the list are doubles and positive numbers, , and if limit is not crossed
        return is_number_array && !overflow && !underflow && is_all_double && is_each_element_positive && !any_element_cross_limit;;
    }
    else {
        print_debug("'char *field_mode' in 'is_field_list_valid()' function is invalid. Please check the code and enter one of the following: 'int_list', 'int_list+', 'double_list', 'double_list+'.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
} 

bool is_there_any_errors(field param, char *additional_error_messages) {
    bool errors = false;
    // Message if parameter is missing
    if (param.missing == true && param.optional == false) {
        printf(COLOR_YELLOW "=> " COLOR_RED "Field " COLOR_MAGENTA "'%s' " COLOR_RED "not found in file contents.\n" RESET_STYLE, param.name);
        errors = true;
    }
    if (param.duplicate == true) {
        printf(COLOR_YELLOW "=> " COLOR_RED "Field " COLOR_MAGENTA "'%s' " COLOR_RED "is declared more than once in file contents.\n" RESET_STYLE, param.name);
        errors = true;
    }
    if (param.missing == false && param.duplicate == false && param.valid == false) {
        printf(COLOR_YELLOW "=> " COLOR_RED "Field " COLOR_MAGENTA "'%s' " COLOR_RED "declared in file contents is invalid:\n" RESET_STYLE, param.name);
        errors = true;
    }
    if (additional_error_messages != NULL) {
        print_error("%s", additional_error_messages);
    }
    
    return errors;   
}

void continue_or_break_execution(field_group *groups, size_t num_groups) {
    // Iterate over all groups
    for (size_t i = 0; i < num_groups; i++) {
        // Iterate over all fields in the current group
        for (size_t j = 0; j < groups[i].len; j++) {
            // Check the .errors member
            if (groups[i].fields[j].errors) {
                print_exit_prog();
                exit(EXIT_FAILURE);
            }
        }
    }
    return;
}

/* ============= FUNCTIONS TO GET INFORMATION FROM SYSTEM ============== */
size_t get_number_of_system_parameters(char *func_name, ode_system *system_list, size_t system_list_len) {
    size_t npar = 0;
    // Check if func_name is NULL
    if (func_name != NULL) {
        // Loop through the system_list to find func_name
        for (int i = 0; i < system_list_len; i++) {
            // check if the system_list[i].func_name is equal to func_name
            if (strcmp(system_list[i].func_name, func_name) == 0) {
                // If true, return the npar (mumber of system parameters) variable
                npar = system_list[i].npar;
            }
        }
    } else {
        print_debug("char *func_name in 'get_number_of_system_parameters()' function NULL. Check the code or input file.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Check also if func_name is not found
    if (npar == 0) {
        print_debug("char *func_name in 'get_number_of_system_parameters()' could not be found in system_list struct. Check the code or input file.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }

    return npar;
}

size_t get_dimension_of_system(char *func_name, ode_system *system_list, size_t system_list_len) {
    size_t dim = 0;
    // Check if func_name is NULL
    if (func_name != NULL) {
        // Loop through the system_list to find func_name
        for (int i = 0; i < system_list_len; i++) {
            // check if the system_list[i].func_name is equal to func_name
            if (strcmp(system_list[i].func_name, func_name) == 0) {
                // If true, return the dim (dimension) variable
                dim = system_list[i].dim;
            }
        }
    } else {
        print_debug("char *func_name in 'get_dimension_of_system()' function NULL. Check the code or input file.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Check also if func_name is not found
    if (dim == 0) {
        print_debug("char *func_name in 'get_dimension_of_system()' could not be found in system_list struct. Check the code or input file.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }

    return dim;
}

/* ==================== FULL CHECKS FOR PARAMETERS ===================== */

void check_fields(char *file_contents, size_t num_fields, field params[num_fields], char **field_names, char *field_mode, size_t *limit_value, bool optional) {
    //print_field_mode(field_mode);
    unsigned int n = 0;
    char *field_content = NULL;
    char *field_value = NULL;
    // Loop through all fields
    for (int i = 0; i < num_fields; i++) {
        // Declare a pointer to memory to hold error messages
        char *error_messages = NULL;
        // Copy field name to params[i].name member 
        if (strlen(field_names[i] + 1) <= INPUT_BUFSIZE) {
            strcpy(params[i].name, field_names[i]);
            params[i].optional = optional;
        }  
        else {
            print_debug("Field name '%s' with length greater than INPUT_BUFSIZE. Reduce the name.\n", field_names[i]);
        }
        // Count field occurrence
        n = count_field_occurrence(file_contents, field_names[i]);
        // Check for missing or duplicate field
        params[i].missing = is_field_missing(n);
        params[i].duplicate = is_field_duplicate(n);
        // Get the corresponding field content based on parameter (field) name and the field value
        field_content = get_field_content(file_contents, params[i].name, field_mode);
        // Check if the field_value is valid based on field_mode
        field_value = get_field_value(field_content);           
        // Check conditions for optional fields
        if (params[i].optional == true && params[i].missing == true) {
            params[i].valid = true;
        }
        else {
            params[i].valid = is_field_value_valid(field_value, limit_value, field_mode, &error_messages); 
        }
        // Check if all conditions are met to assign value to param struct
        if (params[i].missing == false && params[i].duplicate == false && params[i].valid == true) {
            assign_field_value(&params[i], field_value, field_mode);
        }
        // Check for errors
        params[i].errors = is_there_any_errors(params[i], error_messages);
        
        print_field_struct(params[i], field_mode);
        //print_blue("field_content = |%s| | field_value = |%s|\n", field_content, field_value);
        
        // Free Memory
        free(error_messages);
    }
}

void check_list_fields(char *file_contents, size_t num_fields, field params[num_fields], char **field_names, size_t expected_n_cells, bool allow_duplicate_cells, size_t *cell_limit_value, char *field_mode, bool optional) {
    //print_field_mode(field_mode);
    unsigned int n = 0;
    char *field_content = NULL;
    char *field_value = NULL;
    // Loop through all fields
    for (int i = 0; i < num_fields; i++) {
        // Declare a pointer to memory to hold error messages
        char *error_messages = NULL;
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
        // Get the corresponding field content based on parameter (field) name and the field value
        field_content = get_field_content(file_contents, params[i].name, field_mode);
        // Get field value based on field content
        field_value = get_field_list_value(field_content);
        // Get the number of cells of the list    
        size_t n_cells = get_number_of_list_cells_declared(field_value);
        // Store each element of the list into an array of strings
        char **field_values_array = get_each_element_of_list(field_value, n_cells, params[i].name);
        // Check conditions for optional fields
        if (optional == true && params[i].missing == true) {
            params[i].valid = true;
        } 
        else {
            params[i].valid = is_field_list_valid(field_values_array, n_cells, expected_n_cells, allow_duplicate_cells, cell_limit_value, &error_messages, field_mode);
        }
        // Check if all conditions are met to assign values to param struct
        if (params[i].missing == false && params[i].duplicate == false && params[i].valid == true) {
            assign_field_list_values(params, field_values_array, n_cells, field_mode);
        }
        // Check for errors
        params[i].errors = is_there_any_errors(params[i], error_messages);
        print_field_list_struct(params[i], n_cells, field_mode);
        //print_blue("field_content = |%s| | field_value = |%s|\n", field_content, field_value);
    }
} 

/* ============================ MAIN FUNCTION ========================== */

int main(void) {
    /* ============================================================== */
    /* Define Systems and its additional information */
    /* ============================================================== */
    int lorenz_angle_indexes[] = {0, 1};
    ode_system systems[] = {
        // func_name, dim, npar, full_name, group, func(), customcalc_func(), angles, angle_indexes
        {"duffing", 2, 5, "Duffing Oscillator", "HOS", duffing, customcalc_duffing, 0, NULL},
        {"bistable_EH", 3, 8, "Polynomial Bistable Energy Harvester", "HOS", bistable_EH, customcalc_bistable_EH, 0, NULL},                   
        
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
    char *integrator_field_names[] = {"nP", "nDiv", "maxper"};
    size_t integrator_field_names_len = sizeof(integrator_field_names) / sizeof(integrator_field_names[0]);
    char *transient_field_names[] = {"trans"};
    size_t transient_field_names_len = sizeof(transient_field_names) / sizeof(transient_field_names[0]);
    char *IC_field_names[] = {"IC"};
    size_t IC_field_names_len = sizeof(IC_field_names) / sizeof(IC_field_names[0]);
    char *time_field_names[] = {"t0"};
    size_t time_field_names_len = sizeof(time_field_names) / sizeof(time_field_names[0]);
    char *n_RMScalc_field_names[] = {"nrms"};
    size_t n_RMScalc_field_names_len = sizeof(n_RMScalc_field_names) / sizeof(n_RMScalc_field_names[0]);
    char *RMScalc_field_names[] = {"RMS_calc"};
    size_t RMScalc_field_names_len = sizeof(RMScalc_field_names) / sizeof(RMScalc_field_names[0]);
    /* ============================================================== */
    /* Initialize Structs */
    /* ============================================================== */
    field *system_fields = init_field_struct(system_field_names_len, system_field_names);
    field *simulation_fields = init_field_struct(simulation_field_names_len, simulation_field_names);
    field *integrator_fields = init_field_struct(integrator_field_names_len, integrator_field_names);
    field *transient_fields = init_field_struct(transient_field_names_len, transient_field_names);
    field *IC_fields = init_field_struct(IC_field_names_len, IC_field_names);
    field *time_fields = init_field_struct(time_field_names_len, time_field_names);
    field *n_RMScalc_fields = init_field_struct(n_RMScalc_field_names_len, n_RMScalc_field_names);
    field *RMScalc_fields = init_field_struct(RMScalc_field_names_len, RMScalc_field_names);
    /* ============================================================== */
    /* Safety Checks in File Contents */
    /* ============================================================== */
    // Check for system_field_names
    printf(COLOR_YELLOW "=> " RESET_STYLE "Identifying system declared in file contents...\n" RESET_STYLE);
    check_fields(fcontents, system_field_names_len, system_fields, system_field_names, "string", NULL, false);
    // Find system function and customcalc function
    void (*odesys)(int, double *, double, double *, double *) = NULL;
    if (system_fields[1].value != NULL) {
        void (*odesys)(int, double *, double, double *, double *) = find_ode_system(system_fields[1].value, systems, systems_len);
        void (*customcalc)(double *, double *, double, double *, double *, double *, double *, double, int, int, double, int, char **, double *, int) = find_customcalc_func(system_fields[1].value, systems, systems_len);
    }
    // Assign all information of the system to key variables
    size_t npar = get_number_of_system_parameters(system_fields[1].value, systems, systems_len);
    size_t dim = get_dimension_of_system(system_fields[1].value, systems, systems_len);
    // Check for simulation_field_names
    printf(COLOR_YELLOW "=> " RESET_STYLE "Identifying simulation declared in file contents...\n" RESET_STYLE);
    check_fields(fcontents, simulation_field_names_len, simulation_fields, simulation_field_names, "string", NULL, false);
    // Check for remaining fields
    printf(COLOR_YELLOW "=> " RESET_STYLE "Identifying remaining parameters declared in file contents...\n" RESET_STYLE);
    check_fields(fcontents, integrator_field_names_len, integrator_fields, integrator_field_names, "int+", NULL, false);
    if (&((size_t*)integrator_fields->value)[0] != NULL) {
       check_fields(fcontents, transient_field_names_len, transient_fields, transient_field_names, "int+", &((size_t*)integrator_fields->value)[0], false);    
    }   
    check_list_fields(fcontents, IC_field_names_len, IC_fields, IC_field_names, dim, true, NULL, "double_list", false);
    check_fields(fcontents, time_field_names_len, time_fields, time_field_names, "double+", NULL, false);
    // Check for optional fields
    check_fields(fcontents, n_RMScalc_field_names_len, n_RMScalc_fields, n_RMScalc_field_names, "int+", &dim, true);
    if (n_RMScalc_fields->value != NULL && *(size_t*)n_RMScalc_fields->value > 0) {
        size_t limit = dim - 1;
        check_list_fields(fcontents, RMScalc_field_names_len, RMScalc_fields, RMScalc_field_names, *(size_t*)n_RMScalc_fields->value, false, &limit, "int_list+", false);
    }
    // Check for parameter fields 
    char **syspar_field_names = get_names_of_table(npar, "p");
    field *syspar = init_field_struct(npar, syspar_field_names);
    check_fields(fcontents, npar, syspar, syspar_field_names, "double", NULL, false);
    /* ============================================================== */
    /* Check for any errors */
    /* ============================================================== */
    // Create an array containing all field structures
    field_group groups[] = {
        {system_fields, system_field_names_len},
        {simulation_fields, simulation_field_names_len},
        {integrator_fields, integrator_field_names_len},
        {transient_fields, transient_field_names_len},
        {IC_fields, IC_field_names_len},
        {time_fields, time_field_names_len},
        {n_RMScalc_fields, n_RMScalc_field_names_len},
        {RMScalc_fields, RMScalc_field_names_len},
        {syspar, npar}
    };
    // Get the number of field_groups
    size_t num_groups = sizeof(groups) / sizeof(groups[0]);
    // Check if the execution needs to be interrupted based on the validity of the input parameters
    continue_or_break_execution(groups, num_groups);
    /* ============================================================== */
    /* Test System determination */
    /* ============================================================== */
    if (system_fields[1].value != NULL) {
        void (*odesys)(int, double *, double, double *, double *) = find_ode_system(system_fields[1].value, systems, systems_len);
        if (odesys != NULL) {
            num_integrator(odesys, dim, NULL, 0, NULL, NULL);
        } else {
            print_exit_prog();
            exit(EXIT_FAILURE);
        }
    }
    /* ============================================================== */
    /* Free Memory */
    /* ============================================================== */
    free(fcontents);
    free_2D_mem((void**)syspar_field_names, npar);
    free_field_struct(system_field_names_len, system_fields);
    free_field_struct(simulation_field_names_len, simulation_fields); 
    free_field_struct(integrator_field_names_len, integrator_fields); 
    free_field_struct(transient_field_names_len, transient_fields);
    free_field_struct(IC_field_names_len, IC_fields);
    free_field_struct(time_field_names_len, time_fields);
    free_field_struct(n_RMScalc_field_names_len, n_RMScalc_fields);
    free_field_struct(RMScalc_field_names_len, RMScalc_fields);
    free_field_struct(npar, syspar); 
}