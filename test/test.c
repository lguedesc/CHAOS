#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

int main() {

    double *x = malloc(sizeof(*x));
    double *y = malloc(sizeof(y));

    printf("sizeof(*x) = %zu\n", sizeof(*x));
    printf("sizeof(y) = %zu\n", sizeof(y));
    
    char *c = malloc(sizeof(c));
    char *k = malloc(sizeof(k));

    printf("sizeof(*c) = %zu\n", sizeof(*c));
    printf("sizeof(k) = %zu\n", sizeof(k));
    
    return 0;
}
