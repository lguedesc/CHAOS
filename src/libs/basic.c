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

void free_mem(void* first, ...)
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

bool check_if_string_is_number(const char *str, const char *type) {
    int i = 0;
    bool has_decimal = false;
    // Check for optional sign
    if (str[i] == '+' || str[i] == '-')
        i++;
    // Check for digits
    while (str[i] != '\0') {
        // Skip newline and tab characters
        if (str[i] == '\n' || str[i] == '\t') {
            i++;
        }
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