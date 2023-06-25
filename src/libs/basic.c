#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include "msg.h"

int count_int_digits(int number) {
    int count = 0;
    // Handle negative numbers
    if (number < 0) {
        count++;  // Account for the negative sign
        number = abs(number);  // Convert negative number to positive
    }
    // Count the digits
    if (number == 0) {
        count = 1;  // Special case for number 0
    } else {
        while (number != 0) {
            count++;
            number /= 10;
        }
    }
    
    return count;
}

void free_mem(void *first, ...)
{
    va_list args;
    va_start(args, first);
    void* current = first;
    while (current != NULL) {
        free(current);
        current = va_arg(args, void*);
    }
    va_end(args);
}

void free_2D_mem(void **array_2D, int nrows) {
    // Safety check
    if (array_2D == NULL) {
        return;
    }
    // Free memory for each string in the array
    for (int i = 0; i < nrows; i++) {
        free(array_2D[i]);
    }
    // Free memory for the array itself
    free(array_2D);
}

bool check_if_string_is_number_old(const char* str) {
    int i = 0;
    // Check for optional sign
    if (str[i] == '+' || str[i] == '-')
        i++;
    // Check for digits
    while (str[i] != '\0') {
        if (!isdigit(str[i]))
            return false;  // Not a number
        i++;
    }
    return true;  // All characters are numbers
}

bool check_if_string_is_number(const char *str, const char *type, bool only_positive) {
    int i = 0;
    bool has_decimal = false;
    // Check for optional sign
    if (only_positive == true) {
        if (str[i] == '+') {
            i++;
        }
    }
    else {
        if (str[i] == '+' || str[i] == '-') {
            i++;
        }
    }
    // Check for digits
    while (str[i] != '\0') {
        if (isdigit(str[i])) {
            i++;
        }
        else if (str[i] == '.') {
            // Check if decimal point already encountered
            if (has_decimal == true) {
                // More than one decimal point, thus not a number
                return false;   
            }
            has_decimal = true;
            i++;
        }
        else if ((str[i] == '\n') || (str[i] == '\t')) {
            // Allow newline and tab characters
            i++;
        }
        else {
            // Invalid character, not a number
            return false;
        }
    }

    // Check the specifier type
    if (strcmp(type, "int") == 0) {
        // Check if it has a decimal point
        if (has_decimal == true || i == 0) {
            return false;
        }
    }
    else if (strcmp(type, "float") != 0 && (strcmp (type, "double") != 0)) {
        // Invalid type specifier
        print_debug("Invalid type specified in check_if_string_is_number() funtion.\n");
        exit(EXIT_FAILURE);
    }
    return true;
}

bool check_if_string_is_negative_number(const char *str) {
    // Trim leading whitespace
    while (isspace(*str)) {
        str++;
    }
    if (str[0] == '-') {  // Check if the first character is a minus sign
        // Check if the remaining characters are digits
        for (int i = 1; str[i] != '\0'; i++) {
            if (str[i] < '0' || str[i] > '9') {
                return false;  // Found a non-digit character
            }
        }
        return true;  // All characters are digits after the minus sign
    }
    return false;  // First character is not a minus sign
}

char **malloc_string_array(int nstrings, int maxlen) {
    // Create pointer and allocate memory for the number of strings
    char **str_array = malloc(nstrings * sizeof(char*));
    // Safety check
    if (str_array == NULL) {
        print_warning("Failed to allocate memory for the string array\n");
        return NULL;
    }
    // Allocate memory for each string in the array
    for (int i = 0; i < nstrings; i++) {
        str_array[i] = malloc((maxlen + 1) * sizeof(char));
        // Safety check
        if (str_array[i] == NULL) {
            print_warning("Failed to allocate memory for a string in the string array.\n");
            // Clean up the previously allocated memory
            for (int j = 0; j < i; j++) {
                free(str_array[j]);
            }
            free(str_array);
            return NULL;
        }
        // Initialize the string with null characters
        memset(str_array[i], '\0', (maxlen + 1));
    }
    return str_array;
}

void file_safety_check(FILE *file) {
    if (file == NULL) {
        print_error("File to check parameters could not be openeed.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

void ptr_safety_check(void* ptr) {
    if (ptr == NULL) {
        print_debug("NULL pointer found.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}