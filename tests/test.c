#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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

bool is_string_double(char *str) {
    /* THIS FUNCTION SHOULD ONLY BE USED IF *str IS A NUMBER */
    char *ptr;
    // Check if string is not NULL
    if (str != NULL) {
        // Try to convert string to double
        double number = strtod(str, &ptr);
        // Check if is a decimal point in the string
        if (strchr(str, '.') != NULL) {
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

void checkNumberType(char* str) {
    char *ptr;
    double number = strtod(str, &ptr); // try to convert string to double

    // Check if the entire string was converted
    if (ptr == str || *ptr != '\n' && *ptr != '\0') {
        printf("This is not a number.\n");
    } else {
        if (strchr(str, '.') != NULL) { // check if there is a decimal point in the string
            printf("This is a double.\n");
        } else {
            printf("This is an integer.\n");
        }
    }
}

int main() {
    
    char *str = "-10.0";
    char *ptr;

    bool is_number = is_string_number(str);
    printf("is_number = %s\n", is_number ? "true" : "false");
    if (is_number == true) {
        bool is_double = is_string_double(str);
        printf("is_double = %s\n", is_double ? "true" : "false");
    }
    return 0;
}


