#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double getOrder_candidate(double n1, double n2) {
    double diff = fabs(n1 - n2);
    double power;
    if (diff == 0) {
        power = 0;
    }
    else {
        power = log10(diff);
    }
    double order = pow(10, power);
    return order;
}

double getOrder(double n1, double n2) {
    double diff = fabs(n1 - n2);
    double order;
    if (diff == 0) {
        order = 1;
    }
    else {
        order = diff;
    }
    return order;
}

double *get_system_tol_old(int dim, double *xmin, double *xmax) {
    double diff;
    double power;
    double *systol = malloc(dim * sizeof(*systol));
    for (int i = 0; i < dim; i++) {
        diff = fabs(xmax[i] - xmin[i]);
        if (diff == 0) {
            power = 0;
        }
        else {
            power = log10(diff);            
        }        
        systol[i] = pow(10, power);
    }

    return systol;
}

double *get_system_tol(int dim, double *xmin, double *xmax) {
    double diff;
    double *systol = malloc(dim * sizeof(*systol));
    for (int i = 0; i < dim; i++) {
        diff = fabs(xmax[i] - xmin[i]);
        if (diff == 0) {
            systol[i] = 1;
        }
        else {
            systol[i] = diff;            
        }        
    }

    return systol;
}

int main() {
    /*
    double n1, n2;
    
    printf("Enter the first number: ");
    scanf("%lf", &n1);
    
    printf("Enter the second number: ");
    scanf("%lf", &n2);
    
    double order = getOrder(n1, n2);
    printf("%lf - %lf = %lf\n", n1, n2, n1 - n2);
    printf("The order of the set of numbers is %e = %lf\n", order, order);
    */
    int dim = 3;
    double xmax[] = { 0.1, 0, 100 };
    double xmin[] = { -0.32, 0, -50};

    double *systol = get_system_tol(dim, xmin, xmax);

    for(int i = 0; i < dim; i++) {
        printf("systol[%d] = %e = %lf\n", i, systol[i], systol[i]);
    }

    return 0;
}
