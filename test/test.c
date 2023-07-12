#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double value;
    char *name;
    int index;
} cvar;

typedef struct {
    double *x;
    char *par;
    double t0;
    double t;
    double *xrms;
    double *xmin;
    double *xmax;
    double *IC;
    int N;
    int currenttimestep;
    double steadystateperc;
    int ncustomvalues;
    char **customnames;
    size_t maxstrlen;
    double *cvalue;
} cargs;

/*double *customcalc(cvar *var, cargs args) {

}*/

#define getName(var) #var

cvar assign_vars(double variable) {
    cvar newvar;
    // Assign value
    newvar.value = variable;
    // Assign name
    char *variable_name = getName(variable);
    size_t size = sizeof(getName(variable));
    printf("size = %zu\n", size);
    newvar.name = malloc(size + 1);
    strcpy(newvar.name, variable_name);
    printf("newvar.name = %s\n", newvar.name);

    return newvar;
}


int main() {
    
    double Pout = 10*10;
    //double tflip = 2*10; 

    cvar var = assign_vars(Pout);

    printf("var.value = %lf\n", var.value);
    printf("var.name = %s\n", var.name);

    return 0;
}
