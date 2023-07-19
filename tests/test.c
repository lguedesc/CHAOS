#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#define PI (4 * atan(1))
#define TWOPI 2*PI

int main()
{
    
    double xmin = PI*1.1;
    double xmax = 2*PI;

    double xmin_remainder = remainder(xmin, TWOPI);
    double xmax_remainder = remainder(xmax, TWOPI);

    printf("xmin_remainder = %lf\n", xmin_remainder);
    printf("xmax_remainder = %lf\n", xmax_remainder);

    
    return 0;
}
