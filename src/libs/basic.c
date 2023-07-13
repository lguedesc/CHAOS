#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
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

void **alloc_2D_array(int rows, int columns, size_t elementSize) {
    void** array = malloc(rows * sizeof(void*));
    if (array == NULL) {
        print_debug("Memory allocation failed for rows in 'allocate_2D_array()' function.\n");
        return NULL;
    }
    
    for (int i = 0; i < rows; i++) {
        array[i] = malloc(columns * elementSize);
        if (array[i] == NULL) {
            print_debug("Memory allocation failed for columns in 'allocate_2D_array()' function.\n");
            // Free the memory allocated for rows before returning NULL
            for (int j = 0; j < i; j++) {
                free(array[j]);
            }
            free(array);
            return NULL;
        }
    }
    
    return array;
}

char **alloc_string_array(int nstrings, size_t maxlen) {
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

void close_files(int num_files, ...)
{
    va_list args;
    va_start(args, num_files);
    for (int i = 0; i < num_files; i++) {
        FILE* file = va_arg(args, FILE*);
        if (file != NULL) {
            fclose(file);
        }
    }
    va_end(args);
}

void file_safety_check(FILE *file) {
    if (file == NULL) {
        print_error("No such file or directory.\n");
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

void ptr_safety_check(void *ptr, char *ptr_name) {
    if (ptr == NULL) {
        if (ptr_name == NULL) {
            print_debug("'%s' is a NULL pointer, please check the code\n");
        } 
        else {
            print_debug("NULL pointer found, please check the code\n");
        }
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

void double_ptr_safety_check(void** ptr, char *ptr_name) {
    if (ptr == NULL || *ptr == NULL) {
        if (ptr_name == NULL) {
            print_debug("'%s' is a NULL pointer, please check the code\n");
        } 
        else {
            print_debug("NULL pointer found, please check the code\n");
        }
        print_exit_prog();
        exit(EXIT_FAILURE);
    }
}

int get_largest_element_int_array(int *arr, int size) {
    // Check pointer
    ptr_safety_check(arr, "*arr in get_largest_element_int_array()");
    // Assume the first element is the largest
    int largest = arr[0];  
    // Loop through arr to check if any is larger than "largest"
    for (int i = 1; i < size; i++) {
        // If it is larger, update the value of "largest"
        if (arr[i] > largest) {
            largest = arr[i];
        }
    }
    
    return largest;
}