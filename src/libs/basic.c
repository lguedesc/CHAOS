#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

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