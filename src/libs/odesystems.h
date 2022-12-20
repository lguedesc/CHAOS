// Methods
static void error(void);

// General Nonlinear Systems
void lorenz(int dim, double *x, double t, double *par, double *f);
void lotka_volterra_predator_prey(int dim, double *x, double t, double *par, double *f);
void halvorsen(int dim, double *x, double t, double *par, double *f);

// Nonlinear Oscillators
void falksma(int dim, double *x, double t, double *par, double *f);
void vanderpol(int dim, double *x, double t, double *par, double *f);
void pendulum(int dim, double *x, double t, double *par, double *f);
void duffing(int dim, double *x, double t, double *par, double *f);
void linear_oscillator(int dim, double *x, double t, double *par, double *f);
void duffing_2DoF(int dim, double *x, double t, double *par, double *f);
void duffing_vanderpol(int dim, double *x, double t, double *par, double *f);
void linear_oscillator_2DoF(int dim, double *x, double t, double *par, double *f);

// Mechanical Energy Harvesters
void bistable_EH(int dim, double *x, double t, double *par, double *f);
void tristable_EH(int dim, double *x, double t, double *par, double *f);
void pend_oscillator_EH(int dim, double *x, double t, double *par, double *f);
void pend_oscillator_wout_pend_EH(int dim, double *x, double t, double *par, double *f);
void duffing_2DoF_EH(int dim, double *x, double t, double *par, double *f);

// Not Implemented
void duffing_cldyn(int dim, double *x, double t, double *par, double *f);