/*// Methods
static void assign_names(char **strings, const int nvalues, char **names, size_t maxstrlen);
static void spins(double initial_angle, double *previous_angle, double *current_angle, double angle, double *positive_spin, double *negative_spin, int index);
static void time_to_flip(double t, double initial_angle, double current_angle, double *tflip);

// Methods for pend oscillator
static double pend_oscillator_XCM(double X, double rho, double l, double Phi);
static double pend_oscillator_ZCM(double Z, double rho, double l, double Phi);
static double pend_oscillator_dXCM(double dX, double rho, double l, double Phi, double dPhi);
static double pend_oscillator_dZCM(double dZ, double rho, double l, double Phi, double dPhi);
static double pend_oscillator_ddXCM(double ddX, double rho, double l, double Phi, double dPhi, double ddPhi);
static double pend_oscillator_ddZCM(double ddZ, double rho, double l, double Phi, double dPhi, double ddPhi);

// Methods for duffing_2DoF EH
static double duffing_2DoF_EH_XCM(double X1, double X2, double rho);
static double duffing_2DoF_EH_dXCM(double dX1, double dX2, double rho);
static double duffing_2DoF_EH_ddXCM(double ddX1, double ddX2, double rho);
*/
// Custom Calculations
void customcalc(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_bistable_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_pend_oscillator_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
void customcalc_duffing_2DoF_EH(double *x, double *par, double t, double *xrms, double *xmin, double *xmax, double *IC, double t0, int N, int currenttimestep, double steadystateperc, int ncustomvalues, char **customnames, size_t maxstrlen, double *customvalue, int mode);
