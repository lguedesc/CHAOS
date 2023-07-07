
void rk4(int dim, double *x, double t, double h, double *par, double *f, void (*edosys)(int, double *, double, double *, double *)) {
    // Define arrays for temporary storage of the intermediary values of x[dim]
    double tmp[dim], k1[dim], k2[dim], k3[dim], k4[dim];
    // Calculate first slope
    edosys(dim, x, t, par, k1);
    // Assign next value for tmp[dim] to be inserted in the system of edos
    for (int i = 0; i < dim; i++) {
       tmp[i] = x[i] + 0.5 * (h * k1[i]);
    }
    // Calculate second slope
    edosys(dim, tmp, t + 0.5 * h, par, k2);
    // Assign next value for tmp[dim] to be inserted in the system of edos
    for (int i = 0; i < dim; i++) {
        tmp[i] = x[i] + 0.5 * (h * k2[i]);
    }
    // Calculate third slope
    edosys(dim, tmp, t + 0.5 * h, par, k3);
    // Assign next value for tmp[dim] to be inserted in the system of edos
    for (int i = 0; i < dim; i++) {
        tmp[i] = x[i] + h * k3[i];
    }
    // Calculate the fourth slope
    edosys(dim, tmp, t + h, par, k4);
    // Calculate the next value of x[dim]
    for (int i = 0; i < dim; i++) {
        x[i] = x[i] + (h/6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}