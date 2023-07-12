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

double *customcalc(cvar *var, cargs args) {

}



int main() {


    
    return 0;
}
