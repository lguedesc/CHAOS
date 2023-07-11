#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
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