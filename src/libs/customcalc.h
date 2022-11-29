// Methods
static void assign_names(char **strings, const int nvalues, char **names, size_t maxstrlen);
static void spins(double *previous_angle, double *current_angle, double angle, double *positive_spin, double *negative_spin, int index);
static void time_to_flip(double t, double *initial_angle, double current_angle, double *tflip, int *mark, int index);

// Custom Calculations
void customcalc(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_pend_oscillator_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);