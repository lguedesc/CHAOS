#include <stdio.h>
#include <stdbool.h>
#include "src/libs/basic.h"

char **define_IC_names(int dim) {
    
    int digits = 3;
    char **ICnames = malloc_string_array(dim, digits + 3 + 1);


    for (int i = 0; i < dim; i++) {
        sprintf(ICnames[i], "x[%d]", i);
    }
    for (int i = 0; i < dim; i++) {
        printf("ICnames[%d] = %s\n", i, ICnames[i]);
    }

    return ICnames;
}

int main(void) {
    
    int dim = 5;
    char **icnames = define_IC_names(dim);

    free_2D_mem((void**)icnames, dim);
    return 0;
}