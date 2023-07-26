#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#define PI (4 * atan(1))
#define TWOPI 2*PI

int main()
{
    
    double L = 0.005; 
    double g = 9.81;

    double wphi = sqrt(g/L);
    printf("wn = %lf Hz\n", wphi);

    double wz = 39;

    printf("Omega_phi = %lf\n", wphi/wz);
    
    return 0;
}
