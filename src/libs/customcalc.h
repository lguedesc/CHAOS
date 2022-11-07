// Methods
void assign_names(char **strings, const int nvalues, char **names, size_t maxstrlen);

// Custom Calculations
void customcalc(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);