#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

typedef struct {
    int n_angles;
    int *index;
} ang_info;

ang_info *init_angle_struct(int nangles) {
    ang_info *info = malloc(sizeof(ang_info));
    info->n_angles = nangles;
    info->index = malloc(nangles * sizeof(*(info->index)));

    return info;
}

void free_ang_info_struct(ang_info *strct) {
    free(strct->index);
    free(strct);
}

void func(int nangles, ...)
{
    // Declare and initialize variadic argument list
    va_list args;
    va_start(args, nangles);
    // Create ang_info struct with the appropriate size
    ang_info *angle = init_angle_struct(nangles);
    // Assign values to the angle.index array
    for (int i = 0; i < nangles; i++) {
        angle->index[i] = va_arg(args, int);
    }
    va_end(args);
    // Print the array
    for (int i = 0; i < nangles; i++) {
        printf("%d ", angle->index[i]);
    }
    printf("\n");

    // Free dynamically allocated memory
    free_ang_info_struct(angle);
}

int main()
{
    func(5, 1, 2, 3, 4, 5);  // Example usage with angles = 5
    func(3, 6, 7, 8);         // Example usage with angles = 3
    
    double x = 0.0;

    double x_rem = x/6.2342432;
    printf("x_rem = %lf\n", x_rem);


    return 0;
}
