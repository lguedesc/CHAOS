#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI (4 * atan(1))
#define TWOPI 2*PI

int main() {
    
    double xmin4 = -3.3126;
    double xmax4 = 0.411;
    double x4poinc = -1.31;

    double xmin4_rem = remainder(xmin4, TWOPI);
    double xmax4_rem = remainder(xmax4, TWOPI);
    double x4poinc_rem = remainder(x4poinc, TWOPI);

    printf("xmin4 = %lf | remainder = %lf\n", xmin4, xmin4_rem);
    printf("xmax4 = %lf | remainder = %lf\n", xmax4, xmax4_rem);
    printf("x4pnc = %lf | remainder = %lf\n", x4poinc, x4poinc_rem);

    return 0;
}
