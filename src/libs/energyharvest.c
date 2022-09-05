#include <stdio.h>
#include <stdlib.h>
#include "edosystems.h"
#include "iofiles.h"
#include "nldyn.h"

double RMS(double *cum, double measure, int N, int mode) {
    if (mode == 0) {
        // accumulate the value of the square of the measure 
        (*cum) = (*cum) + (measure * measure);
        return (*cum);
    }
    else if (mode == 1) {
        double RMS;
        // Computes the RMS value (cummulative / N times accumulated)
        RMS = sqrt( (*cum) / N );
        return RMS;
    }
    else {
        printf("Failed to compute dimensionless RMS using mode (%d)\n", mode);
        return;
    }
}