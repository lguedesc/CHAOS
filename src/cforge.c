#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <time.h>
#include "libs/msg.h"
#include "libs/basic.h"
#include "libs/iofiles.h"


#define MAX_LINE_LENGTH 10000       // 10 Kb

typedef struct {
    int dim;
    char *name;
    char *outfile;
    char *abrev;
    char *group;
    char *comments;
    char *equations;
    char *lin_equations;
} sys;

sys init_sys_struct(char *name, char *outfile, char *abrev, char *group, char *comments, char *equations, char *lin_equations) {
    sys system;
    // Allocate memory for strings (reduce it later, too big)
    system.name = malloc(strlen(name) + 1);
    system.outfile = malloc(strlen(outfile) + 1);
    system.abrev = malloc(strlen(abrev) + 1);
    system.group = malloc(strlen(group) + 1);
    system.comments = malloc(strlen(comments) + 1);
    system.equations = malloc(strlen(equations) + 1);
    system.lin_equations = malloc(strlen(lin_equations) + 1);
    // Initialize
    system.dim = 0;

    return system;
}

/* File manipulation */

size_t get_file_size(FILE *fp) {
    // Seek the end of the file
    fseek(fp, 0L, SEEK_END);
    // Return the current file position
    size_t size = ftell(fp);
    // Go back to the beggining of the file
    rewind(fp);

    return size;
}

FILE *open_file(char *filename, const char *mode, bool msg) {
    if (msg == true) {
        if (strcmp(mode, "w") == 0) {
        printf("\nCreating new '%s' file...\n", filename);
        }
        else {
            printf("\nSearching for '%s' file...\n", filename);
        }
    }
    FILE *inputfile = fopen(filename, mode);
    if (inputfile == NULL) {
        printf("Could not open file '%s'\n", filename);
        exit(EXIT_FAILURE);
    } 
    else {
        size_t size = get_file_size(inputfile);
        if (msg == true) {
            print_success("File '%s' (%zu bytes) loaded successfully!\n", filename, size);
        }
        return inputfile;
    }
}

size_t get_size_of_longest_line(FILE *file) {
    // Go to the begginning of the file
    rewind(file);
    // Declare variables
    int linenum = 0;
    int current_linenum = 0;
    size_t largest_linesize = 0;
    int c;
    int pos = 0;
    // Reads each character of the file
    while((c = fgetc(file)) != EOF) {
        // increase position
        pos++;
        // check if c char reaches a newline
        if (c == '\n') {
            // increase line number
            current_linenum++;
            if (pos > largest_linesize) {
                largest_linesize = pos;
                linenum = current_linenum;
            }
            // reset position for new line
            pos = 0;
        }
    }
    
    if (largest_linesize > MAX_LINE_LENGTH) {
        print_error("Contents written at line '%d' within the input file exceeds the maximum line length. Please split the information across more lines before running the program again.\n", linenum);
        print_exit_prog();
        exit(EXIT_FAILURE);
    }

    return largest_linesize;
}

/* Check information */
bool continue_program() {
    char choice;
    bool proceed;
    do {
        print_warning("Do you want to continue? [y/n]: ");
        scanf(" %c", &choice);
        if (choice == 'y' || choice == 'Y') {
            // Continue with the program
            proceed = true;
            break;
        } else if (choice == 'n' || choice == 'N') {
            // Exit the program
            proceed = false;
            break;
        } else {
            printf("Invalid choice. Please enter 'y' or 'n'.\n");
        }
    } while (1);
    return proceed;
}

void search_for_key(FILE *file, char *keyword) {
    // Allocate memory for buffer
    size_t buffersize = get_size_of_longest_line(file);
    // Go to the begginning of the file
    rewind(file);
    char *buffer = malloc(buffersize);
    bool found = false;
    int n_found = 0;
    // Seartch for the keyword
    while(fgets(buffer, buffersize, file) != NULL) {
        if (strstr(buffer, keyword) != NULL) {
            found = true;
            n_found++;
        }
    }
    // If keyword not found
    if (found == false) {
        print_error("'%s' key not found within the input file. This information is needed to successfully insert a new system into CHAOS. Please review the input file before running the program again.\n", keyword);
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    if (n_found > 1) {
        print_error("'%s' key found more than once within the input file. Please revise the input file and remove the duplicates before running the program again.\n", keyword);
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // free memory
    free(buffer);
}

void check_section_flags_within_file(FILE *file) {
    char *flags[] = {"# INFO", "# COMMENTS", "# EQUATIONS", "# LINEARIZED EQUATIONS"};
    // get the size of flags[]
    const size_t flags_length = sizeof(flags)/sizeof(flags[0]);
    for (int i = 0; i < flags_length; i++) {
        search_for_key(file, flags[i]);
    }
}

//bool check_if_string_is_number(const char* str) {
//    int i = 0;
    // Check for optional sign
//    if (str[i] == '+' || str[i] == '-')
//        i++;
    // Check for digits
//    while (str[i] != '\0') {
//        if (!isdigit(str[i]))
//            return false;  // Not a number
//        i++;
//    }
//    return true;  // All characters are numbers
//}

void check_group(char *str) {
    if((strcmp(str, "GNL") != 0) && (strcmp(str, "OS") != 0)) {
        print_error("'%s' is a invalid group! Please choose between the options below:\n", str);
        print_error("   1. General Nonlinear Systems by entering 'group = GNL' in the input file;\n");
        print_error("   2. Harmonic Oscillators by entering 'group = OS' in the input file.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

void check_outname(char *str) {
    bool invalid = false;
    char *C_reserved_words[] = {"auto", "break", "case", "char", "const", "continue", "default", "do",
                                  "double", "else", "enum", "extern", "float", "for", "goto", "if", "int",
                                  "long", "register", "return", "short", "signed", "sizeof", "static", 
                                  "struct", "switch", "typedef", "union", "unsigned", "void", "volatile",
                                  "while", "pragma"};
    char *C_first_valid_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_";
    char *C_subsequent_valid_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789";
    // Get the length of the string
    size_t len = strlen(str);
    // Get the length of the C_reserved_words array
    size_t len_reserved = sizeof(C_reserved_words)/sizeof(C_reserved_words[0]); 
    // Check if the name is a reserved keyword
    for (int i = 0; i < len_reserved; i++) {
        if (strcmp(str, C_reserved_words[i]) == 0) {
            invalid = true;
            break;
        }
    }
    // Check if the first character of the string is valid
    size_t C_first_valid_len = strlen(C_first_valid_chars);
    bool first_equal = false;
    for (int i = 0; i < C_first_valid_len; i++) {
        // Check if first character is equal to one of valid chars
        if (str[0] == C_first_valid_chars[i]) {
            first_equal = true;
            // If true, break
            break;
        }
    }
    // Check if the subsequent characters of the string are valid
    size_t C_subsequent_valid_len = strlen(C_subsequent_valid_chars);
    bool subsequent_equal = false;
    // Loop through the string
    for (int i = 1; i < len; i++) {
        // If the iteration is greater than the first check and checker is false after checking a char
        // then it is invalid... break operation
        if ((i > 1) && (subsequent_equal == false)) {
            break;
        }
        // reset checker
        subsequent_equal = false;
        // Loop through the valid chars
        for (int j = 0; j < C_subsequent_valid_len; j++) {
            // Check if the i char of the string is equal to one of the valid chars
            if (str[i] == C_subsequent_valid_chars[j]) {
                subsequent_equal = true;
            }
        }
    }
    // Check if one of the checks are false
    if ((first_equal == false) || (subsequent_equal == false)) {
        invalid = true;
    }

    // determine if the program can continue
    if (invalid == true) {
        print_error("Invalid 'outname' of '%s'. Valid names need to follow the rules listed below:\n", str);
        print_error("   1. The first character of 'outname' is restricted to 'a -> z', 'A -> Z' or '_' (underscore);\n");
        print_error("   2. The subsequent characters of 'outname' are restricted to 'a -> z', 'A -> Z', '0 -> 9' and '_' (underscore);\n");
        print_error("Please, enter a valid 'outname' in the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

int check_dim(char *str) {
    bool valid_dim = check_if_string_is_number(str, "int", true);
    int dim;
    if (valid_dim == true) {
        dim = atoi(str);
        if (dim < 1) {
            print_error("'dim = %d'. 'dim' field should be equal or greater than 1. Please enter a valid 'dim' in the input file before running the program again.\n", dim);
            print_exit_prog();
            exit(EXIT_FAILURE);
        }
    
    } else {
        print_error("'dim' field is not a number. Please enter a valid 'dim' in the input file before running the program again.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    return dim;
}

void check_eqs(char *str, int dim) {
    size_t max_dim_digits = 5;
    // Check if dim is greater than the limit
    if (dim > pow(10, max_dim_digits - 1)) {
        print_error("'dim = %d' field is too large. Please reduce the size of the system or add it manually in the source code.\n", dim);
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Allocate memory for string f[dim]
    size_t func_size = 3 + max_dim_digits;
    char *func = malloc(func_size + 1);
    bool not_found = false;
    bool should_not_be_added = false;
    // Create strings and check for missing fields  
    for (int i = 0; i < dim; i++) {
        snprintf(func, func_size, "f[%d]", i);
        if (strstr(str, func) == NULL) {
            not_found = true;
            print_error("System's dimension 'dim = %d' declared, but equation '%s' could not be found within the input file.\n", dim, func);
        }
    }
    if (not_found == true) {
        print_error("A set of equations within the range f[0] -> f[dim - 1] needs to be added to be compatible with the declared dimension 'dim'.\nInsert the above missing equations below the '# EQUATIONS' field within the input file to successfully insert the new system into CHAOS.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }

    // Checking for possible equations that should not be added
    for (int i = dim; i < max_dim_digits; i++) {
        snprintf(func, func_size, "f[%d]", i);
        if (strstr(str, func) != NULL) {
            should_not_be_added = true;
            print_error("System's dimension 'dim = %d' declared, but equation '%s' was found within the input file.\n", dim, func);
        }
    } 
    if (should_not_be_added == true) {
        print_error("A set of equations within the range f[0] -> f[dim - 1] needs to be added to be compatible with the declared dimension 'dim'.\nRemove the above equations below the '# EQUATIONS' field within the input file to successfully insert the new system into CHAOS.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // free memory
    free(func);
}

void check_lin_eqs(char *str, int dim) {
    size_t max_dim_digits = 5;
    // Check if dim is greater than the limit
    if (dim > pow(10, max_dim_digits - 1)) {
        print_error("'dim = %d' field is too large. Please reduce the size of the system or add it manually in the source code.\n", dim);
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Allocate memory for string f[dim]
    size_t func_size = 7 + max_dim_digits;
    char *func_format_1 = malloc(func_size + 1);
    char *func_format_2 = malloc(func_size + 1);
    bool not_found = false;
    bool should_not_be_added = false;
    // Create strings and check for missing fields  
    for (int i = dim; i < dim + dim*dim; i = i + dim) {
        //print_warning("i = %d\n", i);
        snprintf(func_format_1, func_size, "f[%d + i]", i);
        snprintf(func_format_2, func_size, "f[%d+i]", i);
        if ((strstr(str, func_format_1) == NULL) && (strstr(str, func_format_2) == NULL)) {
            not_found = true;
            print_error("System's dimension 'dim = %d' declared, but equation '%s' could not be found within the input file.\n", dim, func_format_1);
        }
    }
    if (not_found == true) {
        print_error("A set of linearized equations within the range f[dim*1 + i] -> f[dim*dim + i] needs to be added to be compatible with the declared dimension 'dim'.\nInsert the above missing equations below the '# EQUATIONS' field within the input file to successfully insert the new system into CHAOS.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    
    // Checking for possible equations that should not be added
    int num_in_eqs[dim];
    for (int i = 0; i < dim; i++) {
        num_in_eqs[i] = dim*(i+1);
        //print_warning("num_in_eqs = %d\n", num_in_eqs[i]);
    }
    bool equal = false;
    for (int i = 0; i < max_dim_digits + max_dim_digits*max_dim_digits; i++) {
        snprintf(func_format_1, func_size, "f[%d + i]", i);
        snprintf(func_format_2, func_size, "f[%d+i]", i);
        // Check if the number of equation is the number that should be added
        for (int j = 0; j < dim; j++) {
            if (i == num_in_eqs[j]) {
                equal = true;
            }
        }
        if (equal == false) {
            if ((strstr(str, func_format_1) != NULL) || (strstr(str, func_format_2) != NULL)) {
                should_not_be_added = true;
                print_error("System's dimension 'dim = %d' declared, but equation '%s' was found within the input file.\n", dim, func_format_1);
            }
        }
        else {
            equal = false;
        }
    } 
    if (should_not_be_added == true) {
        print_error("A set of equations within the range f[dim*1 + i] -> f[dim*dim + i] needs to be added to be compatible with the declared dimension 'dim'.\nRemove the above equations below the '# LINEARIZED EQUATIONS' field within the input file to successfully insert the new system into CHAOS.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    
    // free memory
    free_mem(func_format_1, func_format_2, NULL);
}

bool check_if_exp_exists(FILE *file, char *exp) {
    bool exists = false;
    // Get the longest line
    size_t buffersize = get_size_of_longest_line(file);
    // Go to the begginning of the file
    rewind(file);
    // Allocate memory for the buffer
    char *buffer = malloc(buffersize);
    while(fgets(buffer, buffersize, file) != NULL) {
        if(strstr(buffer, exp) != NULL) {
            exists = true;
        }
    }
    // Free memory
    free(buffer);
    return exists;
}

void check_conflicts_in_file(char *filepath, char *name, ...) {
    // Open file in read mode
    FILE *file = open_file(filepath, "r", false);
    // Create a variable that flags if conflict exists
    bool conflict = false;
    // Create list of arguments
    va_list args;
    // Declare char pointer to store each argument
    char *str;
    // Initialize va_list object with the first argument
    va_start(args, name);
    // Assign the first argument to the str variable
    str = name;
    // Loop until the value NULL is encountered
    while (str != NULL) {
        bool exists = check_if_exp_exists(file, str);
        if (exists == true) {
            print_error("'%s' declaration already exists in '%s' file. Please choose another name.\n", str, filepath);
            conflict = true;
        }
        // Retrieve the next string argument
        str = va_arg(args, char*);
    }
    // Clean up the va_list object after using it
    va_end(args);
    if (conflict == true) {
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Close file
    fclose(file);
}

void check_info(sys *system, char *name, char *dim, char *outfile, char *abrev, char *group, char *comments, char *equations, char *lin_equations) {
    printf("Checking for conflicting information within the input file...\n");
    // Check if group is valid
    check_group(group);
    // Check if outname is valid for C language
    check_outname(outfile);
    // Initialize struct
    (*system) = init_sys_struct(name, outfile, abrev, group, comments, equations, lin_equations);
    // Check if dim is valid and store remaining information in system struct    
    system->dim = check_dim(dim);
    // Check if equations and lin_equations are valid based on the declared dimension ('dim')
    check_eqs(equations, system->dim);
    check_lin_eqs(lin_equations, system->dim);
    // Assign checked information to struct
    strcpy(system->name, name);
    strcpy(system->outfile, outfile);
    strcpy(system->abrev, abrev);
    strcpy(system->group, group);
    strcpy(system->comments, comments);
    strcpy(system->equations, equations);
    strcpy(system->lin_equations, lin_equations);
    print_success("Success in the initial checking of the input file!\n");
    // Search for conflicts in main.c, odesystems.c, odesystems.h, customcalc.c, customcalc.h, such as equal outnames, etc...
    printf("Checking for conflicts between the input information and CHAOS source code...\n");
    check_conflicts_in_file("src/main.c", system->name, system->outfile, NULL);
    check_conflicts_in_file("src/libs/odesystems.c", system->outfile, NULL);
    check_conflicts_in_file("src/libs/odesystems.h", system->outfile, NULL);
    check_conflicts_in_file("src/libs/customcalc.c", system->outfile, NULL);
    check_conflicts_in_file("src/libs/customcalc.h", system->outfile, NULL);
    print_success("No conflicts found!\n");

    printf("\nINFO:\n");
    printf("name: %s\n", name);
    printf("dim: %s\n", dim);
    printf("outfile: %s\n", outfile);
    printf("abrev: %s\n", abrev);
    printf("group: %s\n", group);
    printf("\n");
    printf("COMMENTS:\n%s\n\n", comments);
    printf("EQUATIONS:\n%s\n\n", equations);
    printf("LINEARIZED EQUATIONS:\n%s\n\n", lin_equations);
    print_warning("The above information will be added to CHAOS.\n");
    // check to continue
    bool proceed = continue_program();
    if (proceed == false) {
        print_exit_prog();
        exit(EXIT_SUCCESS);
    }
}

/* String Manipulation */
void remove_unwanted_whitespace(char *str) {
    char *start = str;
    char *end = str + strlen(str) - 1;
    // Trim leading whitespace
    while (isspace(*start)) {
        start++;
    }
    // Trim trailing whitespace
    while (end > start && isspace(*end)) {
        end--;
    }
    // Null-terminate the trimmed string
    *(end + 1) = '\0';
    // Shift the trimmed string to the beginning of the input string
    memmove(str, start, end - start + 2);
}

char *add_identation(char *str, int n) {
    int len = strlen(str);
    int i, j;
    int newline_count = 0;
    // Count the number of newline characters in the input string
    for (i = 0; i < len; i++) {
        if (str[i] == '\n') {
            newline_count++;
        }
    }
    // Calculate the new length of the modified string with tabs added (tabs in the beggining of each line)
    int new_len = len + n*(newline_count + 1);  // Each newline will be replaced by '\n\t'
    // Allocate memory for the modified string
    char *modified_str = malloc((new_len + 1) * sizeof(char));
    if (modified_str == NULL || n <= 0) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }
     // Insert tab at the beginning of the string
    int k;
    for (k = 0; k < n; k++) {
        modified_str[k] = '\t';
    }
    for (i = 0, j = k; i < len; i++, j++) {
        if (str[i] == '\n') {
            // Add the '\n' to the j modified_string position and add 1 to j (j=j+1)
            modified_str[j++] = '\n';
            // Add tab to the position next to '\n'
            if (n == 1) {
                // If n = 1, need only 1 operation
                modified_str[j] = '\t';
            } 
            else {
                // If n > 1
                for (int m = 0; m < n; m++) {
                    // Check if it is the last tab
                    if (m == n - 1)
                        modified_str[j] = '\t';
                    else {
                        // If not the last tab, put a tab and advance one j position
                        modified_str[j++] = '\t';
                    }            
                }        
            }
        } else {
            modified_str[j] = str[i];   
        }
    }   
    // Add null character at the end of the modified string
    modified_str[new_len] = '\0';  

    return modified_str;
}

/* Store information */
char *store_string_block(FILE *file, char *keyword) {
    size_t buffersize = get_size_of_longest_line(file);    
    // Allocate memory for a buffer
    char *buffer = malloc(buffersize);
    // Declare control variables
    bool key_found = false;
    int init_block_pos = 0;
    int final_block_pos = 0;
    char *string_block;
    int pos = 0;
    bool hash_found = false;
    // Go to the beggining of the file 
    rewind(file);
    // Start reading line by line and storing into buffer
    while(fgets(buffer, buffersize, file) != NULL) {
        // Get the current file position
        pos = ftell(file);
        if (key_found == false) {
            // Identify position in which content block begins
            if (strstr(buffer, keyword) != NULL) {
                init_block_pos = pos;
                key_found = true;
            }
        }
        /* Check if the keyword was found and if the position of the stream is different
           from the keyword position */
        if ((key_found == true) && (pos > init_block_pos)) {
            // Search for the "#" keyword and 
            if (strstr(buffer, "#")) {
                // determine position of beggining of the line that contains '#'
                final_block_pos = pos - strlen(buffer);
                hash_found = true;
                break;
            }
        }
    }
    // Go the end of the file if key is found and '#' is not found, and get final position
    if ((key_found == true) && (hash_found == false)) {
        fseek(file, 0, SEEK_END);
        final_block_pos = ftell(file);
    }

    if (key_found == false) {
        print_error("Could not find '%s' field within the input file. Please add the field before continuing.", keyword);
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // free memory
    free(buffer);
    // Move to the initial block position
    if (fseek(file, init_block_pos, SEEK_SET) != 0) {
        print_error("DEBUG ERROR: Error setting file position.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Calculate the size of the data to be read
    long block_size = final_block_pos - init_block_pos;
    // Insert the string in a memory block   
    if (block_size != 0) {
        // Allocate memory for the block size
        string_block = malloc(block_size + 1);
        // Read data from file
        size_t bytes_read = fread(string_block, 1, block_size, file);
        // Check if any bytes were read
        if (bytes_read > 0) {
            // Add null terminator at the end of the string
            string_block[bytes_read] = '\0';
        }
    }
    else {
        // Allocate memory for the block size
        string_block = malloc(2);
        strcpy(string_block, " \0");
    }

    return string_block;
}

char *get_assigned_value(FILE *file, char *input_string, char *key) {
    size_t buffersize = get_size_of_longest_line(file);
    char *delim = NULL; 
    char *line = malloc(buffersize);
    char *value = NULL;
    // Search for the first occurrence of key in the input_string and creater a pointer to it
    char *buffer = strstr(input_string, key);
    // Check if key was declared in the input file
    if (buffer != NULL) {
        // get line that key is defined
        for (int i = 0; i < strlen(buffer); i++) {
            line[i] = buffer[i]; 
            if(buffer[i] == '\n') {
                line[i] = '\0';
                break;
            }            
        }
        if (line != NULL) {
            // If line is found, search for the delimiters '=' or ':'
            delim = strchr(line, '='); 
            if (delim == NULL) {
                delim = strchr(line, ':');
            }
            // If the delimiter is found
            if (delim != NULL) {
                // get the value after the delimiter
                value = delim + 1;
                // remove whitespaces
                remove_unwanted_whitespace(value);
            }
        }
        
    }
    else {
        print_error("Missing '%s' parameter in the input file. Please insert '%s' before running the program again.\n", key, key);        
        print_exit_prog();        
        exit(EXIT_FAILURE);
    }

    return value;
}

/* Add information */
void add_ode(sys system) {
    char *filepath_c = "src/libs/odesystems.c";
    // Open library file to read and append system of equations
    FILE *cfile = open_file(filepath_c, "ab+", true);
    printf("Writing contents to '%s' file...\n", filepath_c);
    // Write contents .c file
    fprintf(cfile, "\nvoid %s(int dim, double *x, double t, double *par, double *f) {\n", system.outfile);
    fprintf(cfile, "\t/*\n");
    fprintf(cfile, "%s\n", system.comments);
    fprintf(cfile, "\t*/\n");
    fprintf(cfile, "\tif (dim == %d) {\n", system.dim);
    fprintf(cfile, "%s\n", system.equations);
    fprintf(cfile, "\t}\n");
    fprintf(cfile, "\telse if (dim == %d) {\n", system.dim*system.dim + system.dim);   
    fprintf(cfile, "%s\n", system.equations);
    fprintf(cfile, "\t\tfor(int i = 0; i < %d; i++) {\n", system.dim);
    fprintf(cfile, "%s\n", system.lin_equations);
    fprintf(cfile, "\t\t}\n");
    fprintf(cfile, "\t}\n");
    fprintf(cfile, "\telse {\n");
    fprintf(cfile, "\t\terror();\n");
    fprintf(cfile, "\t}\n");
    fprintf(cfile, "}\n");
    // Close file
    fclose(cfile);
    // Message
    print_success("Contents written successfully in '%s'!\n", filepath_c);

    char *filepath_h = "src/libs/odesystems.h";
    // Open header file to read and append the function declaration of system of equations
    FILE *hfile = open_file(filepath_h, "ab+", true);
    printf("Writing contents to '%s' file...\n", filepath_h);
    // Write content
    fprintf(hfile, "\nvoid %s(int dim, double *x, double t, double *par, double *f);\n", system.outfile);
    // Close file
    fclose(hfile);
    // Message
    print_success("Contents written successfully in '%s'!\n", filepath_h);

}

void add_customcalc(sys system) {
    // Open library c file to read and append system of equations
    char *filepath_c = "src/libs/customcalc.c";
    FILE * cfile = open_file(filepath_c, "ab+", true);
    // Write content in c file
    printf("Writing contents to '%s' file...\n", filepath_c); 
    fprintf(cfile, "\nvoid customcalc_%s(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode) {\n", system.outfile);
    fprintf(cfile, "\treturn;\n");
    fprintf(cfile, "}\n");
    // Close file c file
    fclose(cfile);
    // message
    print_success("Contents written successfully in '%s'!\n", filepath_c);

    // Open header file to read and append the function declaration of system of equations
    char *filepath_h = "src/libs/customcalc.h";
    FILE *hfile = open_file(filepath_h, "ab+", true);
    // Write content
    printf("Writing contents to '%s' file...\n", filepath_h);
    fprintf(hfile, "\nvoid customcalc_%s(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);\n", system.outfile);
    // Close file
    fclose(hfile);
    print_success("Contents written successfully in '%s'!\n", filepath_h);

}

void add_number_of_systems(char* group) {
    char *filepath = "src/main.c";
    // Open main file to read and write the information
    FILE *file = open_file(filepath, "rb+", true);
    // Allocate memory for the buffer and for the keyword that will be searched for within main file
    size_t buffersize = get_size_of_longest_line(file);
    char *buffer = malloc(buffersize);
    char *keyword = malloc(buffersize);
    printf("Updating number of '%s' to '%s' file...\n", group, filepath);
    // Search for keyword to be changed 
    sprintf(keyword, "#define NUM_OF_%s_SYSTEMS", group);
    int number = 0;
    long current_position = 0;
    bool found = false;
    // Go to the beggining of the file
    rewind(file);
    // Search for the keyword and get the current number of systems declared in main file
    while (fgets(buffer, buffersize, file) != NULL) {
        //printf("buffer = %s\n", buffer);
        if (strstr(buffer, keyword) != NULL) {
            found = true;
            sscanf(buffer, "%*s %*s %d\n", &number);
            // Get also the current position of the line that NUM_OF_SYSTEMS is declared
            current_position = ftell(file) - strlen(buffer);
            break;
        }
    }
    // If it doesnt find the keyword, give error and exit program.
    if (found == false) {
        char *error = malloc(MAX_LINE_LENGTH * sizeof(char));
        sprintf(error,"'%s' expression not found in the main.c file.\nThe main.c file can be corrupted.\n", keyword);
        reset_program(error);
        free(error);
    }
    // Increase the number by one
    number++;
    // Move to the position of the line
    fseek(file, current_position, SEEK_SET);
    // Replace the line with the updated number
    fprintf(file, "#define NUM_OF_%s_SYSTEMS %d", group, number);
    // free memory
    free_mem(keyword, buffer, NULL);
    // close file
    print_success("File updated successfuly!\n");
    printf("Closing file...\n");
    fclose(file);
}

int occurrence_of_expr(FILE *file, char *expr, long *last_position, long *line_position) {
    // Variable that holds the number of a expression occurs in a file
    int n = 0;
    bool found = false;
    // Allocate memory for buffer
    size_t buffersize = get_size_of_longest_line(file);
    char *buffer = malloc(buffersize);
    // Go to the begginning of the file
    rewind(file);
    // Search for the keyword and get the current number
    while(fgets(buffer, buffersize, file) != NULL) {
        // If it finds an occurrence, get the position in the file
        if (strstr(buffer, expr) != NULL) {
            found = true;
            n++;
            // Get positions
            (*last_position) = ftell(file);
            (*line_position) = (*last_position) - strlen(buffer);
        }
    }
    // If it doesnt find the search_keyword, give error and exit program.
    if (found == false) {
        char *error = malloc(MAX_LINE_LENGTH * sizeof(char));
        sprintf(error,"'%s' expression not found in the main.c file.\nThe main.c file may be corrupted.\n", expr);
        reset_program(error);
    }    
    // Free memory
    free(buffer);    
    // Return the number of times it encountered the expr
    return n;
}

char *split_and_add_text(FILE *file, long split_pos, char *str_add) {
    // Find the file size
    size_t filesize = get_file_size(file);
    // Go to the beginning of the file 
    rewind(file);
    // Check if split_pos is within valid range
    if (split_pos < 0 || split_pos > filesize) {
        print_error("DEBUG Error: Invalid split position.\n");
        exit(EXIT_FAILURE);
    }
    // Allocate memory for temporary buffers
    size_t remainingsize = filesize - split_pos; 
    char *temp_top = malloc(split_pos + 1);
    char *temp_bot = malloc(remainingsize + 1);
    // Read the content above split_pos
    fseek(file, 0, SEEK_SET);
    fread(temp_top, sizeof(char), split_pos, file);
    temp_top[split_pos] = '\0';
    // Read the content below split_pos
    fseek(file, split_pos, SEEK_SET);
    fread(temp_bot, sizeof(char), remainingsize, file);
    temp_bot[remainingsize] = '\0';
    // Calculate the size of the new string
    size_t newsize = strlen(temp_top) + strlen(str_add) + strlen(temp_bot);
    // Allocate memory for the new string
    char *result = malloc(newsize + 1);
    // Concatenate the strings in the desired order
    strcpy(result, temp_top);
    strcat(result, str_add);
    strcat(result, temp_bot);
    // free memory
    free_mem(temp_top, temp_bot, NULL);
    
    return result;
}

void add_system_information(sys system) {
    char *filepath = "src/main.c";
    // Open main file to read and write the information
    FILE *file = open_file(filepath, "rb+", true);
    size_t maxlinesize = get_size_of_longest_line(file);
    // Go to the begginning of the file
    rewind(file);
    // Defines the word to be searched for in the file
    char *keyword = malloc(maxlinesize);
    sprintf(keyword, "#define %s_OUTPUTNAME_", system.group);
    // Locate the last occurrence of keyword and position
    long last_position = 0;
    long line_position = 0;
    int number = occurrence_of_expr(file, keyword, &last_position, &line_position);
    // add number of the new system
    number++;
    // Construct the new content that will be inserted in the file
    char *newline_1 = malloc(MAX_LINE_LENGTH);
    char *newline_2 = malloc(MAX_LINE_LENGTH);
    char *newline_3 = malloc(MAX_LINE_LENGTH);
    sprintf(newline_1, "\n#define %s_FUNC_%d %s\n", system.group, number, system.outfile);
    sprintf(newline_2, "#define %s_CUSTOM_%d customcalc_%s\n", system.group, number, system.outfile);
    sprintf(newline_3, "#define %s_OUTPUTNAME_%d \"%s\"\n", system.group, number, system.outfile);    
    size_t new_size = strlen(newline_1) + strlen(newline_2)  + strlen(newline_3);
    char *add_string = malloc(new_size + 1); 
    strcpy(add_string, newline_1);
    strncat(add_string, newline_2, strlen(newline_2));
    strncat(add_string, newline_3, strlen(newline_3));
    // Reconstruct information
    char* new_content = split_and_add_text(file, last_position, add_string);
    // Close file
    printf("Closing old %s file...\n", filepath);
    fclose(file);    
    // Create new empty file
    FILE *newfile = open_file(filepath, "wb", true);
    // Write new content to file 
    fwrite(new_content, sizeof(char), strlen(new_content), newfile);    
    // Close file
    fclose(newfile);
    // free memory
    free_mem(new_content, keyword, newline_1, newline_2, newline_3, NULL);
}

void add_system_name(sys system) {
    char *filepath = "src/main.c";
    // Open main file to read and write the information
    FILE *file = open_file(filepath, "rb+", true);
    size_t maxlinesize = get_size_of_longest_line(file);
    // Go to the begginning of the file
    rewind(file);
    // Defines the word to be searched for in the file
    char *keyword = malloc(MAX_LINE_LENGTH);
    sprintf(keyword, "char *%ssystemNames", system.group);
    // Locate the last occurrence of keyword and position
    long last_position = 0;
    long line_position = 0;
    int number = occurrence_of_expr(file, keyword, &last_position, &line_position);    
    // Go to the beggining of the function and find its end
    char *buffer = malloc(maxlinesize);
    fseek(file, last_position, SEEK_SET);
    bool end_found = false;
    char *end = "};";
    while(fgets(buffer, maxlinesize, file) != NULL) {
        if(strstr(buffer, end) != NULL) {
            last_position = ftell(file);
            end_found = true;
            break;
        }
    }
    // If it doesnt find the search_keyword, give error and exit program.
    if (end_found == false) {
        char *error = malloc(MAX_LINE_LENGTH);
        sprintf(error,"'%s' expression not found in the main.c file.\nThe main.c file may be corrupted.\n", end);
        reset_program(error);
    }
    // Inverse loop through the file to find the first occurrence of character "
    char ch;
    for (int i = last_position; i > 0; i--) {
        // Go to the previous 2 characters
        fseek(file, -2, SEEK_CUR);
        // Advance one character and check
        ch = fgetc(file);
        if (ch == '\"') {
            last_position = ftell(file);
            break;
        }
    }
    // Construct the string to be added in file
    char *str1 = ",\n\t\t\t\t\t\t\t\t\t\t\t\"";
    char *str2 = "\"";
    size_t new_strlen = strlen(str1) + strlen(system.name) + strlen(str2); 
    char *add_string = malloc(new_strlen + 1);    
    strcpy(add_string, str1);
    strncat(add_string, system.name, strlen(system.name));
    strncat(add_string, str2, strlen(str2));
    // Split text, add new information and reconstruct text
    char *new_content = split_and_add_text(file, last_position, add_string);
    // Close file
    printf("Closing old %s file...\n", filepath);
    fclose(file);    
    // Create new empty file
    FILE *newfile = open_file(filepath, "wb", true);
    // Write new content to file 
    fwrite(new_content, sizeof(char), strlen(new_content), newfile);    
    // Close file
    fclose(newfile);
    // free memory
    free_mem(keyword, new_content, buffer, add_string, NULL);
}

void add_system_call(sys system) {
    char *filepath = "src/main.c";
    // Open main file to read and write the information
    FILE *file = open_file(filepath, "rb+", true);
    size_t maxlinesize = get_size_of_longest_line(file);
    // Go to the begginning of the file
    rewind(file);
    // Defines the word to be searched for in the file
    char *keyword = malloc(MAX_LINE_LENGTH);
    sprintf(keyword, "execute_%s_modules(module", system.group);
    // Locate the last occurrence of keyword and position
    long last_position = 0;
    long line_position = 0;
    int number = occurrence_of_expr(file, keyword, &last_position, &line_position);
    // Find the next "default" expression in the main file
    char *buffer = malloc(maxlinesize);
    fseek(file, line_position, SEEK_SET);
    bool end_found = false;
    char *end = "default:";
    while(fgets(buffer, maxlinesize, file) != NULL) {
        if(strstr(buffer, end) != NULL) {
            last_position = ftell(file);
            end_found = true;
            break;
        }
    }
    // If it doesnt find the search_keyword, give error and exit program.
    if (end_found == false) {
        char *error = malloc(MAX_LINE_LENGTH);
        sprintf(error,"'%s' expression not found in the main.c file.\nThe main.c file may be corrupted.\n", end);
        reset_program(error);
    }
    // Inverse loop through the file to find the first occurrence of character ';'
    char ch;
    for (int i = last_position; i > 0; i--) {
        // Go to the previous 2 characters
        fseek(file, -2, SEEK_CUR);
        // Advance one character and check
        ch = fgetc(file);
        if (ch == ';') {
            last_position = ftell(file);
            break;
        }
    }
    // Construct the string to be added in file
    number++;                                       // Increase number to add case
    int ndigits = count_int_digits(number);
    char *str1 = "\n\t\tcase ";
    char str2[ndigits + 1];
    sprintf(str2, "%d", number);
    char *str3 = ":\n\t\t\texecute_";
    char *str4 = "_modules(module, ";
    char *str5 = "_FUNC_";
    char *str6 = ", ";
    char *str7 = "_CUSTOM_";
    char *str8 = "_OUTPUTNAME_";
    char *str9 = "systemNames[";
    int number_arr = number - 1;
    int narrdigits = count_int_digits(number_arr);
    char str10[narrdigits + 1];
    sprintf(str10, "%d", number_arr);
    char *str11 = "]);\n\t\t\tbreak;\n";
    // Check for system group to consider the differences
    size_t new_strlen;
    char *add_string = NULL;
    if (strcmp(system.group, "OS") == 0) {
        new_strlen = strlen(str1) + strlen(str2) + strlen(str3) + strlen(system.group) + 
                     strlen(str4) + strlen(system.group) + strlen(str5) + strlen(str2) +
                     strlen(str6) + strlen(system.group) + strlen(str7) + strlen(str2) + 
                     strlen(str6) + strlen(system.group) + strlen(str8) + strlen(str2) + 
                     strlen(str6) + strlen(system.group) + strlen(str9) + strlen(str10) + 
                     strlen(str11);
        // Allocate proportional memory
        add_string = malloc(new_strlen + 1);
        // Construct the string
        strcpy(add_string, str1);
        strncat(add_string, str2, strlen(str2));
        strncat(add_string, str3, strlen(str3));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str4, strlen(str4));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str5, strlen(str5));
        strncat(add_string, str2, strlen(str2));
        strncat(add_string, str6, strlen(str6));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str7, strlen(str7));
        strncat(add_string, str2, strlen(str2));
        strncat(add_string, str6, strlen(str6));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str8, strlen(str8));
        strncat(add_string, str2, strlen(str2));
        strncat(add_string, str6, strlen(str6));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str9, strlen(str9));
        strncat(add_string, str10, strlen(str10));
        strncat(add_string, str11, strlen(str11));
    } 
    else if (strcmp(system.group, "GNL") == 0) {
        new_strlen = strlen(str1) + strlen(str2) + strlen(str3) + strlen(system.group) + 
                     strlen(str4) + strlen(system.group) + strlen(str5) + strlen(str2) + 
                     strlen(str6) + strlen(system.group) + strlen(str8) + strlen(str2) + 
                     strlen(str6) + strlen(system.group) + strlen(str9) + strlen(str10) + 
                     strlen(str11);
        // Allocate proportional memory
        add_string = malloc(new_strlen + 1);
        // Construct the string
        strcpy(add_string, str1);
        strncat(add_string, str2, strlen(str2));
        strncat(add_string, str3, strlen(str3));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str4, strlen(str4));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str5, strlen(str5));
        strncat(add_string, str2, strlen(str2));
        strncat(add_string, str6, strlen(str6));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str8, strlen(str8));
        strncat(add_string, str2, strlen(str2));
        strncat(add_string, str6, strlen(str6));
        strncat(add_string, system.group, strlen(system.group));
        strncat(add_string, str9, strlen(str9));
        strncat(add_string, str10, strlen(str10));
        strncat(add_string, str11, strlen(str11));
    }
    else {
        print_error("'%s' is a invalid group! Please choose between the options below:\n", system.group);
        print_error("   1. General Nonlinear Systems by entering 'group = GNL' in the input file;\n");
        print_error("   2. Harmonic Oscillators by entering 'group = OS' in the input file.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
    // Split text, add new information and reconstruct text
    char *new_content = split_and_add_text(file, last_position, add_string);
    // Close file
    printf("Closing old %s file...\n", filepath);
    fclose(file);
    // Create new empty file
    FILE *newfile = open_file(filepath, "wb", true);
    // Write new content to file 
    fwrite(new_content, sizeof(char), strlen(new_content), newfile);    
    // Close file
    fclose(newfile);
    // free memory
    free_mem(keyword, new_content, buffer, add_string, NULL);
}

/* Welcome and Tutorial */
void file_tutorial() {
    printf("\nBelow there is an example of how a input file must be formatted and organized:\n\n");
    print_purple("# INFO:\n");
    print_purple("dim = 2\n");
    print_purple("name = Test Oscillator\n");
    print_purple("abrev = TOsc\n");
    print_purple("outfile = test_oscillator\n");
    print_purple("group = OS\n");
    print_purple("\n# COMMENTS:\n");
    print_purple("par[0] = Forcing frequency\n");
    print_purple("par[1] = Forcing amplitude\n");
    print_purple("par[2] = Damping Coefficient\n");
    print_purple("par[3] = Spring Stiffness\n");
    print_purple("\n# EQUATIONS:\n");
    print_purple("f[0] = x[1];\n");
    print_purple("f[1] = -(1/par[2])*(par[2]*x[1] + par[3]*x[0]) + par[1]*sin(par[0]*t);\n");
    print_purple("\n# LINEARIZED EQUATIONS:\n");
    print_purple("f[2 + i] = x[4 + i];\n");
    print_purple("f[4 + i] = -2*par[2]*par[3]*x[4 + i] - par[3]*par[3]*x[2 + i];\n\n");

    printf("1. The '# INFO:' block holds the information about how CHAOS will handle\n");
    printf("   your system. It require the following fields:\n");
    printf("\t- dim:     The dimension of the dynamical system, related to\n"); 
    printf("\t           the amount of first order ordinary differential\n"); 
    printf("\t           equations (ODEs) of the system;\n");
    printf("\t- name:    The name of the dynamical system that will be printed\n");
    printf("\t           in the simulation information files generated by the\n");
    printf("\t           CHAOS package;\n");
    printf("\t- abrev:   The abbreviation of the dynamical system name used by\n");
    printf("\t           the CHAOS package to identify simulation requirements;\n");
    printf("\t- outfile: The base name of the output files that will contain\n");
    printf("\t           the simulation results and information generated by\n");
    printf("\t           the CHAOS package;\n");
    printf("\t- group:   The group the custom dynamical system belongs.\n");
    printf("\t           Currently, the CHAOS package supports two types of\n");
    printf("\t           dynamical systems: General Nonlinear Systems (GNL)\n");
    printf("\t           and Harmonic Oscillators (OS). This field is very\n");
    printf("\t           important as each group have its own unique available\n");
    printf("\t           set of features within the CHAOS package. Also, some\n");
    printf("\t           features as bifurcation diagrams and Poincare maps\n");
    printf("\t           are computed differently in each group.\n\n");
    printf("2. The '# COMMENTS:' block can be used as annotations about the system\n");
    printf("   that will be stored into the CHAOS source code. For example, it can\n");
    printf("   hold information about the parameters of the system related to the\n");
    printf("   implementation of the equations. The '# COMMENTS:' block is not \n");
    printf("   mandatory and can be left empty.\n\n");
    printf("3. The '# EQUATIONS:' block holds the implementation of the system of\n");
    printf("   equations in the form of 'f = dx/dt = g(x)', where 'x[]' are the state\n");
    printf("   space variables of the system, 'par[]' are its constant parameters,\n");
    printf("   't' is the time, and 'f[]' is the derivative of the state variables.\n\n");
    printf("4. The '# LINEARIZED EQUATIONS:' block holds the implementation of the\n");
    printf("   linearized equations in the form of 'f = dg/dt = g(x)'. This field is\n");
    printf("   used to calculate the Lyapunov spectra. The method utilized in the CHAOS\n");
    printf("   package utilizes 'dim' sets of linearized equations of the system,\n");
    printf("   resulting in 'n = dim*dim' linearized equations. To simplify, these\n");
    printf("   equations are declared in the form of sets, being 'i = {0, ..., dim}'\n");
    printf("   the number of the set. For that, 'f[]' and 'x[]' are declared adding\n");
    printf("   'i' in the array index. For example, for a system with 'dim = 2',\n");
    printf("   there will be 2 linearized equations within each set declared as f[dim*1 + i]\n");
    printf("   and f[dim*2 + i]. Now, for a system with 'dim = 3', there will be 3\n");
    printf("   linearized equations within each set declared as f[dim*1 + i], f[dim*2 + i]\n");
    printf("   and f[dim*3 + i], and so on. x[] array follows the same rule. To find\n");
    printf("   these linearized equations one must solve the system lf = [A]*lx,\n");
    printf("   where [A] is the Jacobian Matrix of the system, 'lx' is the vector\n");
    printf("   containing all the new state space variables of the linearized equations,\n");
    printf("   and 'lf' if the vector containing all the derivatives of 'lx'vector.\n\n"); 
    printf("   Blocks '# EQUATIONS:' and '# LINEARIZED EQUATIONS:' are mandatory and\n");
    printf("   must follow the C language standards, including the ';' at the end of\n");
    printf("   each equation.\n\n");

}

bool asks_for_tutorial() {
    char choice;
    bool proceed;
    printf("To begin, CHAOS Forge requires an input file with the some key information.\n");
    print_warning("Do you need help to create your input file? [y/n]: ");
    do {
        scanf(" %c", &choice);
        if (choice == 'y' || choice == 'Y') {
            // Continue with the program
            proceed = true;
            break;
        } else if (choice == 'n' || choice == 'N') {
            // Exit the program
            proceed = false;
            break;
        } else {
            printf("Invalid choice. Please enter 'y' or 'n'.\n");
        }
    } while (1);
    return proceed;
}

void welcome() {
    print_blue("\nWelcome to CHAOS Forge! v.1.0\n\n");
    print_blue("This tool allows you to expand the capabilities of the CHAOS package by adding\n");
    print_blue("your own custom dynamical systems, allowing CHAOS to meet your specific needs.\n\n");                
    bool tutorial = false;
    tutorial = asks_for_tutorial();
    if (tutorial == true) {
        file_tutorial();
    }
}

/* Backup */
bool asks_for_backup() {
    char choice;
    bool proceed;
    print_warning("It is highly recommended to make a backup before continuing.\n");
    print_warning("Do you want to create a backup? [y/n]: ");
    do {
        scanf(" %c", &choice);
        if (choice == 'y' || choice == 'Y') {
            // Continue with the program
            proceed = true;
            break;
        } else if (choice == 'n' || choice == 'N') {
            // Exit the program
            proceed = false;
            break;
        } else {
            printf("Invalid choice. Please enter 'y' or 'n'.\n");
        }
    } while (1);
    return proceed;
}

void copy_file(char *filepath, char *dest) {
    size_t bytes_read;
    // Open source file
    FILE *file = open_file(filepath, "rb", false);
    // Open destination file
    FILE *destfile = open_file(dest, "wb", false);
    // Get the size pf
    size_t filesize = get_file_size(file);
    char *buffer = malloc(filesize);    
    // Copy the contents of the source file to the destination file
    while ((bytes_read = fread(buffer, 1, sizeof(buffer), file)) > 0) {
        fwrite(buffer, 1, bytes_read, destfile);
    }
    // free memory
    free(buffer);
    // Close files
    fclose(file); fclose(destfile);
}   

void get_current_date_time(int *day, char month[4], int *year, int *hour, int *minute, int *second) {
    time_t current_time;
    struct tm *local_time;

    int mon;
    // Get the current time
    current_time = time(NULL);

    // Convert the current time to the local time zone
    local_time = localtime(&current_time);

    // Retrieve the date components
    (*day) = local_time->tm_mday;
    mon = local_time->tm_mon + 1;  // Month starts from 0
    (*year) = local_time->tm_year + 1900;  // Year starts from 1900

    // Retrieve the time components
    (*hour) = local_time->tm_hour;
    (*minute) = local_time->tm_min;
    (*second) = local_time->tm_sec;

    // Convert month to string form
    if (mon == 1) { strcpy(month, "jan"); }
    else if (mon == 2) { strcpy(month, "feb"); }
    else if (mon == 3) { strcpy(month, "mar"); }
    else if (mon == 4) { strcpy(month, "apr"); }
    else if (mon == 5) { strcpy(month, "may"); }
    else if (mon == 6) { strcpy(month, "jun"); }
    else if (mon == 7) { strcpy(month, "jul"); }
    else if (mon == 8) { strcpy(month, "aug"); }
    else if (mon == 9) { strcpy(month, "sep"); }
    else if (mon == 10) { strcpy(month, "oct"); }
    else if (mon == 11) { strcpy(month, "nov"); }
    else if (mon == 12) { strcpy(month, "dez"); }
    else {
        print_error("Error retrieving month to create backup.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

char *create_bkp_dir() {
    char *rawdir = "backup/";
    char *dir = convert_dir(rawdir);
    // Get time and date
    int day; char mon[4]; int year; int hour; int min; int sec;
    get_current_date_time(&day, mon, &year, &hour, &min, &sec);
    //Create full string directory
    size_t n = 22;
    char *fulldir = malloc(n);
    snprintf(fulldir, strlen(dir) + n + 1, "%s%d-%s-%d_%dh%dm%ds/", dir, day, mon, year, hour, min, sec);
    // Check if directory exists. If not, create
    directory_exists(fulldir);
    return fulldir;
}

void backup() {
    // Asks for backup
    bool proceed = asks_for_backup();
    if (proceed != true) {
        return;
    }
    // Create directory
    char *backup_dir = create_bkp_dir();
    // Source Directories
    char *main_fulldir = convert_dir("src/main.c");
    char *ode_c_fulldir = convert_dir("src/libs/odesystems.c");
    char *ode_h_fulldir = convert_dir("src/libs/odesystems.h");
    char *custom_c_fulldir = convert_dir("src/libs/customcalc.c");
    char *custom_h_fulldir = convert_dir("src/libs/customcalc.h");
    // Copy files
    char *srcfilenames[5] = {main_fulldir, ode_c_fulldir, ode_h_fulldir, custom_c_fulldir, custom_h_fulldir};
    char *destfilenames[5] = { "main.c", "odesystems.c", "odesystems.h", "customcalc.c", "customcalc.h" };
    size_t size = strlen(backup_dir) + 12 + 1;
    char destdir[size]; 
    for (int i = 0; i < 5; i++) {
        snprintf(destdir, size, "%s%s", backup_dir, destfilenames[i]);
        copy_file(srcfilenames[i], destdir);
    }
}

int main(void) {
    /* 0. Welcome and instructions */
    welcome();
    /* 1. Read input file */
    char *filename = get_input_filename();
    FILE *file = open_file(filename, "rb", true);
    /* 2. Check if input file is missing any flag field or INFO field*/
    printf("Checking if input file '%s' have missing or duplicate any section flag...\n", filename);
    check_section_flags_within_file(file);
    print_success("Section flags successfully declared!\n");
    printf("Checking if input file '%s' have missing or duplicate any INFO field...\n", filename);
    search_for_key(file, "name");
    search_for_key(file, "outfile");
    search_for_key(file, "dim");
    search_for_key(file, "abrev");
    search_for_key(file, "group");
    print_success("INFO fields successfully declared!\n");
    /* 3. Store relevant information of the file in strings */
    sys system;
    // INFO block
    //printf("Reading INFO group.\n");
    char *info_block = store_string_block(file, "INFO");
    char *outfile = get_assigned_value(file, info_block, "outfile");
    char *dim = get_assigned_value(file, info_block, "dim");
    char *name = get_assigned_value(file, info_block, "name");
    char *abrev = get_assigned_value(file, info_block, "abrev");
    char *group = get_assigned_value(file, info_block, "group");
    // COMMENTS block 
    char *comments_block = store_string_block(file, "COMMENTS");
    remove_unwanted_whitespace(comments_block);
    // EQUATIONS block 
    char *equations_block = store_string_block(file, "EQUATIONS");
    remove_unwanted_whitespace(equations_block);
    // LINEARIZED EQUATIONS block
    char *lin_equations_block = store_string_block(file, "LINEARIZED EQUATIONS");
    remove_unwanted_whitespace(lin_equations_block);    
    // Close file
    fclose(file);   
    /* 4. Check if information stored in strings are valid and assign to system struct.
          Also check for any conflicts between input information and source code */
    check_info(&system, name, dim, outfile, abrev, group, comments_block, equations_block, lin_equations_block);
    /* 5. Format information */
    system.comments = add_identation(system.comments, 2);
    system.equations = add_identation(system.equations, 2);
    system.lin_equations = add_identation(system.lin_equations, 3);
    /* 6. Make a backup for files that will be modified */
    backup();
    /* 7. Modify CHAOS source code */
    add_ode(system);
    add_customcalc(system);
    add_number_of_systems(system.group);
    add_system_information(system);
    add_system_name(system);
    add_system_call(system);
    // free memory
    free_mem(filename, info_block, comments_block,
             equations_block, lin_equations_block, NULL);
}