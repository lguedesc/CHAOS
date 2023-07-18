#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>


void init_ptr(double **ptr, size_t rows, size_t cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            ptr[i][j] = 0;
        }
    }
}

void init_x(double *x, size_t cols) {
    for (int j = 0; j < cols; j++) {
        x[j] = 0;
    }
}

void third_level(double **ptr, double *x, size_t rows, size_t cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            ptr[i][j] = 1;
        }
    }
    for (int j = 0; j < cols; j++) {
        x[j] = 1;
    }
}

void second_level(double **ptr, double *x, size_t rows, size_t cols) {
    third_level(ptr, x, rows, cols);
}

void first_level(double **ptr, double *x, size_t rows, size_t cols) {
    second_level(ptr, x, rows, cols);
}

int main()
{
    // Allocate double pointer
    size_t rows = 5, cols = 2;
    double **ptr = malloc(rows * sizeof (*ptr));
    for (int i = 0; i < rows; i++) {
        ptr[i] = malloc(cols * sizeof (**ptr));
    }
    // Allocate pointer
    double *x = malloc(cols * sizeof (*x));
    // Initialize pointers
    init_ptr(ptr, rows, cols);
    init_x(x, cols);
    // Print initialized pointers
    printf("Initialized Pointers:\n\n");
    for (int j = 0; j < cols; j++) {
        printf("x[%d] = %lf\n", j, x[j]);
    }
    printf("\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("ptr[%d][%d] = %lf\n", i, j, ptr[i][j]);
        }
    }
    printf("\n\n");

    // call functions 
    first_level(ptr, x, rows, cols);
    printf("After Function Call Pointers:\n\n");
    for (int j = 0; j < cols; j++) {
        printf("x[%d] = %lf\n", j, x[j]);
    }
    printf("\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("ptr[%d][%d] = %lf\n", i, j, ptr[i][j]);
        }
    }
    printf("\n");
    


    // Free memory
    for (int i = 0; i < rows; i++) {
        free(ptr[i]);
    }
    free(ptr);
    free(x);
    return 0;
}
