// Custom Calculations
void customcalc(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_pend_oscillator_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_duffing_2DoF_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_linear_2DoF_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_adeodato_sma_oscillator(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);

void customcalc_chuas_circuit(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
