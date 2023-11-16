#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_SYS_DECLARATIONS 3

#define HOS_FUNC_1 duffing
#define HOS_FUNC_2 duffing_2DoF
#define HOS_FUNC_3 vanderpol

#define HOS_CUSTOM customcalc

#define HOS_OUTPUTNAME_1 "duffing"
#define HOS_OUTPUTNAME_2 "duffing_2DoF"
#define HOS_OUTPUTNAME_3 "vanderpol"

#define HOS_DIM_1 2
#define HOS_DIM_2 4
#define HOS_DIM_3 2

#define HOS_NPAR_1 5
#define HOS_NPAR_2 9
#define HOS_NPAR_3 3

#define HOS_ANGLES_1 0
#define HOS_ANGLES_2 0
#define HOS_ANGLES_3 0

// Helper macros for stringification and concatenation
#define STR(x) #x
#define STRINGIFY(x) STR(x)
#define CONCAT(x, y) x##y

// Function to compare the input string with HOS_FUNC_N
#define COMPARE_SYSTEM(i) \
    if (strcmp(system, STRINGIFY(HOS_FUNC_##i)) == 0) { \
        printf("\nMatch found for system: %s\n", system); \
        printf("Custom function: %s\n", STRINGIFY(HOS_CUSTOM)); \
        *outname = strdup(STRINGIFY(HOS_OUTPUTNAME_##i)); \
        *dim = CONCAT(HOS_DIM_, i); \
        *npar = CONCAT(HOS_NPAR_, i); \
        printf("Output name: %s\n", *outname); \
        printf("Dimension: %zu\n", *dim); \
        printf("Npar: %zu\n", *npar); \
        printf("Angles: %d\n", CONCAT(HOS_ANGLES_, i)); \
        return; \
    }

void compareSystem(const char *system, size_t *dim, size_t *npar, char **outname) {
    COMPARE_SYSTEM(1);
    COMPARE_SYSTEM(2);
    COMPARE_SYSTEM(3);
    printf("\nNo match found for system: %s\n", system);
}

int main() {
    char *system = "duffing";
    size_t dim, npar;
    char *outname;
    compareSystem(system, &dim, &npar, &outname);
    printf("Assigned values: Dim=%zu, Npar=%zu, Outname=%s\n", dim, npar, outname);

    system = "vanderpol";
    compareSystem(system, &dim, &npar, &outname);
    printf("Assigned values: Dim=%zu, Npar=%zu, Outname=%s\n", dim, npar, outname);

    system = "lorenz"; // Not in the list
    compareSystem(system, &dim, &npar, &outname);
    printf("Assigned values: Dim=%zu, Npar=%zu, Outname=%s\n", dim, npar, outname);

    // Don't forget to free the memory allocated for outname
    free(outname);

    return 0;
}
