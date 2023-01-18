#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/* To be added */

// Double Pendulum
// Falk SMA 2DoF Oscillator
// Nonsmooth Oscillator

/* To be added EH */
// Nonsmooth EH

static void error(void) {
    printf("  ERROR:\n");
    printf("  Wrong dimension allocated for the system of equations.\n");
    printf("  Please check the input file.\n");
    printf("  Exiting Program...\n");
    exit(1);
}

// General Nonlinear Systems
void halvorsen(int dim, double *x, double t, double *par, double *f) {
    // OMEGA = par[0]
    // a = par[1]

    if (dim == 3) {
        f[0] = -par[1]*x[0] - 4.0*x[1] - 4.0*x[2] - x[1]*x[1];
        f[1] = -par[1]*x[1] - 4.0*x[2] - 4.0*x[0] - x[2]*x[2];
        f[2] = -par[1]*x[2] - 4.0*x[0] - 4.0*x[1] - x[0]*x[0];
    } 
    else if (dim == 12) {
        f[0] = -par[1]*x[0] - 4.0*x[1] - 4.0*x[2] - x[1]*x[1];
        f[1] = -par[1]*x[1] - 4.0*x[2] - 4.0*x[0] - x[2]*x[2];
        f[2] = -par[1]*x[2] - 4.0*x[0] - 4.0*x[1] - x[0]*x[0];
        for (int i = 0; i < 3; i++) {   
            f[3 + i] = -par[1]*x[3+i] + (-4.0 - 2.0*x[1])*x[6+i] - 4.0*x[9+i];
            f[6 + i] = -par[1]*x[6+i] + (-4.0 - 2.0*x[2])*x[9+i] - 4.0*x[3+i];
            f[9 + i] = -par[1]*x[9+i] + (-4.0 - 2.0*x[0])*x[3+i] - 4.0*x[6+i];
        }
    }
    else {
        error();
    }    
}

void lorenz(int dim, double *x, double t, double *par, double *f) {
    // OMEGA = par[0]
    // sigma = par[1]
    // rho = par[2]
    // beta = par[3]
    if (dim == 3) {
        f[0] = par[1]*(x[1] - x[0]);
        f[1] = x[0]*(par[2] - x[2]) - x[1];
        f[2] = x[0]*x[1] - par[3]*x[2];    
    } 
    else if (dim == 12) {
        f[0] = par[1]*(x[1] - x[0]);
        f[1] = x[0]*(par[2] - x[2]) - x[1];
        f[2] = x[0]*x[1] - par[3]*x[2];
        for (int i = 0; i < 3; i++) {
            f[3 + i] = par[1]*(x[6 + i] - x[3 + i]);
            f[6 + i] = (par[2] - x[2])*x[3 + i] - x[6 + i] - x[0]*x[9 + i];
            f[9 + i] = x[1]*x[3 + i] + x[0]*x[6 + i] - par[3]*x[9 + i];
        }
    }
    else {
        error();
    }    
}

void lotka_volterra_predator_prey(int dim, double *x, double t, double *par, double *f) {
    // OMEGA   = par[0]   |   x[0] = x => number of prey      
    // alpha   = par[1]   |   x[1] = y => number of predator  
    // beta    = par[2]   |        
    // delta   = par[3]   |      
    // gamma   = par[4]   |           
    
    if (dim == 2) {
        f[0] = par[1]*x[0] - par[2]*x[0]*x[1];
        f[1] = par[3]*x[0]*x[1] - par[4]*x[1]; 
    } 
    else if (dim == 6) {
        f[0] = par[1]*x[0] - par[2]*x[0]*x[1];
        f[1] = par[3]*x[0]*x[1] - par[4]*x[1];
        for (int i = 0; i < 2; i ++) {
            // Add later
            f[2 + i] = (par[1] - par[2]*x[1])*x[2+i] - par[2]*x[0]*x[4+i];
            f[4 + i] = (-par[4] + par[3]*x[0])*x[4+i] + par[3]*x[1]*x[2+i];
        }
    }
    else {
        error();
    }
}

// Nonlinear Oscillators
void falksma(int dim, double *x, double t, double *par, double *f) {
    // OMEGA   = par[0]
    // gamma   = par[1]
    // zeta    = par[2]
    // theta   = par[3]
    // beta    = par[4]
    // theta_A = par[5]
    if (dim == 2) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - par[2]*x[1] - (par[3] - 1)*x[0] + par[4]*x[0]*x[0]*x[0] - (par[4]*par[4] / (4*(par[5] - 1)))*x[0]*x[0]*x[0]*x[0]*x[0];
    }
    else if (dim == 6) {
        f[0] = x[1];
        f[1] = f[1] = par[1]*sin(par[0] * t) - par[2]*x[1] - (par[3] - 1)*x[0] + par[4]*x[0]*x[0]*x[0] - (par[4]*par[4] / (4*(par[5] - 1)))*x[0]*x[0]*x[0]*x[0]*x[0];
        for (int i = 0; i < 2; i++) {
            f[2 + i] = x[4 + i];
            f[4 + i] = -par[2]*x[4 + i] - (1 + 3*par[4]*x[0]*x[0] - par[3] -  (5/4)*(par[4]*par[4]/(4*(1 - par[5])))*x[0]*x[0]*x[0]*x[0])*x[2 + i];
        }
    }
    else {
        error();
    }
}

void vanderpol(int dim, double *x, double t, double *par, double *f) {
    // OMEGA = par[0]
    // gamma = par[1]
    // mu    = par[2]
    if (dim == 2) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) + par[2]*(1 - x[0]*x[0])*x[1] - x[0];
    }
    else if (dim == 6) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) + par[2]*(1 - x[0]*x[0])*x[1] - x[0];
        for (int i = 0; i < 2; i++) {
            f[2 + i] = x[4 + i];
            f[4 + i] = (-1 - 2*par[2]*x[0]*x[1])*x[2 + i] + par[2]*(1 - x[0]*x[0])*x[4 + i];
        }
    }
    else {
        error();
    }
}

void pendulum(int dim, double *x, double t, double *par, double *f) {
    // OMEGA = par[0]
    // gamma = par[1]
    // zeta  = par[2]
    // g     = par[3]
    // L     = par[4]
    if (dim == 2) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - par[2]*x[1] - (par[3]/par[4])*sin(x[0]);
    }
    else if (dim == 6) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - par[2]*x[1] - (par[3]/par[4])*sin(x[0]);
        for (int i = 0; i < 2; i++) {
            f[2 + i] = x[4 + i];
            f[4 + i] = -par[2]*x[4 + i] - (par[3]/par[4])*cos(x[0])*x[2 + i];
        } 
    }
    else {
        error();
    }
}

void duffing(int dim, double *x, double t, double *par, double *f) {
    // OMEGA = par[0]
    // gamma = par[1]
    // zeta = par[2]
    // alpha = par[3]
    // beta = par[4]
    
    if (dim == 2) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0];    
    } 
    else if (dim == 6) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0];
        for (int i = 0; i < 2; i++) {
            f[2 + i] = x[4 + i];
            f[4 + i] = -par[3]*x[2 + i] - 3*par[4]*x[0]*x[0]*x[2 + i] - 2*par[2]*x[4 + i];
        }
    }
    else {
        error();
    }    
}

void linear_oscillator(int dim, double *x, double t, double *par, double *f) {
    // OMEGA = par[0]
    // gamma = par[1]
    // zeta = par[2]
    // omega_n = par[3]
    if (dim == 2) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*par[3]*x[1] - par[3]*par[3]*x[0];
    }
    else if (dim == 6) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*par[3]*x[1] - par[3]*par[3]*x[0];
        for (int i = 0; i < 2; i++) {
            f[2 + i] = x[4 + i];               
            f[4 + i] = -2*par[2]*par[3]*x[4 + i] - par[3]*par[3]*x[2 + i];       
        }
    }
    else {
        error();
    }
}

void duffing_2DoF(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA  = par[0] | alpha2 = par[5]
       gamma  = par[1] | beta1  = par[6]
       zeta1  = par[2] | beta2  = par[7]
       zeta2  = par[3] | rho    = par[8]
       alpha1 = par[4] |              */
    if (dim == 4) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] + 2*par[3]*(x[3] - x[1]) - par[4]*x[0] - par[6]*x[0]*x[0]*x[0] + par[5]*(x[2] - x[0]) + par[7]*pow((x[2] - x[0]), 3);
        f[2] = x[3];
        f[3] = -(1/par[8])*(2*par[3]*(x[3] - x[1]) + par[5]*(x[2] - x[0]) + par[7]*pow((x[2] - x[0]), 3));
    }	
    else if (dim == 20) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] + 2*par[3]*(x[3] - x[1]) - par[4]*x[0] - par[6]*x[0]*x[0]*x[0] + par[5]*(x[2] - x[0]) + par[7]*pow((x[2] - x[0]), 3);
        f[2] = x[3];
        f[3] = -(1/par[8]) * (2*par[3]*(x[3] - x[1]) + par[5]*(x[2] - x[0]) + par[7]*pow((x[2] - x[0]), 3));
        for (int i = 0; i < 4; i++) {
            f[4 + i] = x[8 + i];
            f[8 + i] = (-par[4] - par[5] - 3*par[6]*x[0]*x[0] - 3*par[7]*pow((x[2] - x[0]), 2))*x[4 + i] 
                       + (par[5] + 3*par[7]*pow((x[2] - x[0]), 2))*x[12 + i] + (-2*par[2] - 2*par[3])*x[8 + i] + 2*par[3]*x[16 + i];
            f[12 + i] = x[16 + i];
            f[16 + i] = (1/par[8])*((-par[5] - 3*par[7]*pow((x[2] - x[0]), 2))*x[12 + i] + (par[5] + 3*par[7]*pow((x[2] - x[0]), 2))*x[4 + i]
			            + 2*par[3]*x[8 + i] - 2*par[3]*x[16 + i]);
        }
    }
    else {
        error();
    }
}

void duffing_vanderpol(int dim, double *x, double t, double *par, double *f) {
    // OMEGA   = par[0]   |   zeta = par[5]      
    // gamma   = par[1]   |     
    // epsilon = par[2]   |        
    // alpha   = par[3]   |      
    // beta    = par[4]   |           
    
    if (dim == 2) {
        f[0] = x[1];
        f[1] = par[1]*par[0]*par[0]*sin(par[0]*t) - par[2]*x[1]*(x[0]*x[0] - 1) - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] - par[5]*x[0]*x[0]*x[0]*x[0]*x[0]; 
    } 
    else if (dim == 6) {
        f[0] = x[1];
        f[1] = par[1]*par[0]*par[0]*sin(par[0]*t) - par[2]*x[1]*(x[0]*x[0] - 1) - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] - par[5]*x[0]*x[0]*x[0]*x[0]*x[0]; 
        for (int i = 0; i < 2; i ++) {
            // Add later
            f[2 + i] = x[4 + i];
            f[4 + i] = (-par[3] - 3*par[4]*x[0]*x[0] - 2*par[2]*x[0]*x[1] - 5*par[5]*x[0]*x[0]*x[0]*x[0])*x[2 + i] - par[2]*(x[0]*x[0] - 1)*x[4 + i];
        }
    }
    else {
        error();
    }
}

void linear_oscillator_2DoF(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA   = par[0]   |   OMEGA_s = par[5]   |   x1  = x[0] | 
       gamma   = par[1]   |                      |   dx1 = x[1] |    
       rho     = par[2]   |                      |   x2  = x[2] |
       zeta_1  = par[3]   |                      |   dx2 = x[3] |
       zeta_2  = par[4]   |                      |              |      */
    if (dim == 4) { 
        f[0] = x[1];
        f[1] = par[1]*par[0]*par[0]*sin(par[0]*t) - 2*par[3]*x[1] + 2*par[4]*(x[3] - x[1]) - x[0] + par[2]*par[5]*par[5]*(x[2] - x[0]);
        f[2] = x[3];
        f[3] = par[1]*par[0]*par[0]*sin(par[0]*t) - (2/par[2])*par[4]*(x[3] - x[1]) - par[5]*par[5]*(x[2] - x[0]);
    }
    else if (dim == 42) {
        f[0] = x[1];
        f[1] = par[1]*par[0]*par[0]*sin(par[0]*t) - 2*par[3]*x[1] + 2*par[4]*(x[3] - x[1]) - x[0] + par[2]*par[5]*par[5]*(x[2] - x[0]);
        f[2] = x[3];
        f[3] = par[1]*par[0]*par[0]*sin(par[0]*t) - (2/par[2])*par[4]*(x[3] - x[1]) - par[5]*par[5]*(x[2] - x[0]);
        for (int i = 0; i < 4; i ++) {
            f[4 + i] = 0;
            f[8 + i] = 0;
            f[12 + i] = 0;
            f[16 + i] = 0;
        }
    }
    else {
        error();
    }
}

// Mechanical Energy Harvesters
void bistable_EH(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA = par[0]   |   chi    = par[5]  
       gamma = par[1]   |   varphi = par[6]
       zeta  = par[2]   |   kappa  = par[7]
       alpha = par[3]   |
       beta  = par[4]   |                */
    if (dim == 3) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] + par[5]*x[2];
        f[2] = -par[6]*x[2] - par[7]*x[1];
    }
    else if (dim == 12) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] + par[5]*x[2];
        f[2] = -par[6]*x[2] - par[7]*x[1];
        for (int i = 0; i < 3; i ++) {
            f[3 + i] = x[6 + i];
            f[6 + i] = -par[3]*x[3 + i] - 3*par[4]*x[0]*x[0]*x[3 + i] - 2*par[2]*x[6 + i] + par[5]*x[9 + i];
            f[9 + i] = -par[6]*x[9 + i] - par[7]*x[6 + i];
        }
    }
    else {
        error();
    }
}

void tristable_EH(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA = par[0]   |   sigma  = par[5]  
       gamma = par[1]   |   chi    = par[6]
       zeta  = par[2]   |   varphi = par[7]
       alpha = par[3]   |   kappa  = par[8]
       beta  = par[4]   |                */
    if (dim == 3) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] - par[5]*x[0]*x[0]*x[0]*x[0]*x[0] + par[6]*x[2];
        f[2] = -par[7]*x[2] - par[8]*x[1];
    }
    else if (dim == 12) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0] - par[5]*x[0]*x[0]*x[0]*x[0]*x[0] + par[6]*x[2];
        f[2] = -par[7]*x[2] - par[8]*x[1];
        for (int i = 0; i < 3; i ++) {
            f[3 + i] = x[6 + i];
            f[6 + i] = - 2*par[2]*x[6 + i] + par[6]*x[9 + i] - (par[3] + 3*par[4]*x[0]*x[0] + 5*par[5]*x[0]*x[0]*x[0]*x[0])*x[3 + i];
            f[9 + i] = -par[7]*x[9 + i] - par[8]*x[6 + i];
        }
    }
    else {
        error();
    }
}

void pend_oscillator_EH(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA   = par[0]   |   zeta_z    = par[5]   |   l         = par[10]   |   chi_PZ = par[15]       |   x[0] = x       |   x[5] = dphi/dt
       gamma   = par[1]   |   zeta_t    = par[6]   |   varphi_PZ = par[11]   |   chi_EM = par[16]       |   x[1] = dx/dt   |   x[6] = v
       mu      = par[2]   |   OMEGA_s   = par[7]   |   kappa_PZ  = par[12]   |                          |   x[2] = z       |   x[7] = i
       rho     = par[3]   |   OMEGA_phi = par[8]   |   varphi_EM = par[13]   |                          |   x[3] = dz/dt   |
       zeta_x  = par[4]   |   OMEGA_t   = par[9]   |   kappa_EM  = par[14]   |                          |   x[4] = phi     |                   */       
    const double pi = 4 * atan(1);
    if (dim == 8) { 
        f[0] = x[1];
        f[1] = (1/(1 + par[3]))*(-(1 + par[3]*cos(x[4])*cos(x[4]))*(2*par[4]*x[1] + par[7]*par[7]*x[0]) + (par[3]/2)*(2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*sin(x[4]))
                + (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*cos(x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*sin((pi/180)*par[2]); 
        f[2] = x[3];
        f[3] = (1/(1 + par[3]))*(-(1 + par[3]*sin(x[4])*sin(x[4]))*(2*par[5]*x[3] + x[2] - par[15]*x[6]) + (par[3]/2)*(2*par[4]*x[1] + par[7]*par[7]*x[0])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*cos(x[4]))
                - (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*sin(x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*cos((pi/180)*par[2]);
        f[4] = x[5];
        f[5] = -(1 + par[3])*(2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7]) + (1/par[10])*((2*par[4]*x[1] + par[7]*par[7]*x[0])*cos(x[4]) - (2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(x[4]));
        f[6] = -par[12]*x[3] - par[11]*x[6];
        f[7] = -par[14]*x[5] - par[13]*x[7];
    } 
    else if (dim == 72) {
        f[0] = x[1];
        f[1] = (1/(1 + par[3]))*(-(1 + par[3]*cos(x[4])*cos(x[4]))*(2*par[4]*x[1] + par[7]*par[7]*x[0]) + (par[3]/2)*(2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*sin(x[4]))
                + (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*cos(x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*sin((pi/180)*par[2]); 
        f[2] = x[3];
        f[3] = (1/(1 + par[3]))*(-(1 + par[3]*sin(x[4])*sin(x[4]))*(2*par[5]*x[3] + x[2] - par[15]*x[6]) + (par[3]/2)*(2*par[4]*x[1] + par[7]*par[7]*x[0])*sin(2*x[4]) + par[3]*par[10]*x[5]*x[5]*cos(x[4]))
                - (2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7])*par[3]*par[10]*sin(x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*cos((pi/180)*par[2]);
        f[4] = x[5];
        f[5] = -(1 + par[3])*(2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4]) - par[16]*x[7]) + (1/par[10])*((2*par[4]*x[1] + par[7]*par[7]*x[0])*cos(x[4]) - (2*par[5]*x[3] + x[2] - par[15]*x[6])*sin(x[4]));
        f[6] = -par[12]*x[3] - par[11]*x[6];
        f[7] = -par[14]*x[5] - par[13]*x[7];
        for (int i = 0; i < 8; i ++) {
            f[8 + i] = x[16 + i];
            f[16 + i] = (-(par[7]*par[7]*(1 + par[3]*cos(x[4])*cos(x[4]))*x[8 + i]) - 2*par[4]*(1 + par[3]*cos(x[4])*cos(x[4]))*x[16 + i] + par[3]*cos(x[4])*sin(x[4])*x[24 + i] + par[5]*par[3]*sin(2*x[4])*x[32 + i] + 
                         par[3]*(sin(2*x[4])*(par[7]*par[7]*x[0] + 2*par[4]*x[1]) + par[10]*cos(x[4])*((1 + par[3])*par[9]*par[9] + x[5]*x[5]) + 
                         cos(2*x[4])*(par[10]*(1 + par[3])*par[8]*par[8] + x[2] + 2*par[5]*x[3] - par[15]*x[6]) - par[10]*(1 + par[3])*sin(x[4])*(par[9]*par[9]*x[4] + 2*par[6]*x[5] - par[16]*x[7]))*x[40 + i] + 
                         2*par[10]*par[3]*(par[6]*(1 + par[3])*cos(x[4]) + sin(x[4])*x[5])*x[48 + i] - par[3]*par[15]*cos(x[4])*sin(x[4])*x[56 + i] - par[10]*par[3]*(1 + par[3])*par[16]*cos(x[4])*x[64 + i])/(1 + par[3]);
            f[24 + i] = x[32 + i];
            f[32 + i] = (par[3]*par[7]*par[7]*cos(x[4])*sin(x[4])*x[8 + i] + par[4]*par[3]*sin(2*x[4])*x[16 + i] - (1 + par[3]*sin(x[4])*sin(x[4]))*x[24 + i] - 2*par[5]*(1 + par[3]*sin(x[4])*sin(x[4]))*x[32 + i] + 
                         par[3]*(cos(2*x[4])*(par[7]*par[7]*x[0] + 2*par[4]*x[1]) - sin(x[4])*(par[10]*((1 + par[3])*par[9]*par[9] + x[5]*x[5]) + 
                         2*cos(x[4])*(par[10]*(1 + par[3])*par[8]*par[8] + x[2] + 2*par[5]*x[3] - par[15]*x[6])) - par[10]*(1 + par[3])*cos(x[4])*(par[9]*par[9]*x[4] + 2*par[6]*x[5] - par[16]*x[7]))*x[40 + i] - 
                         2*par[10]*par[3]*(par[6]*(1 + par[3])*sin(x[4]) - cos(x[4])*x[5])*x[48 + i] + (par[15] + par[3]*par[15]*sin(x[4])*sin(x[4]))*x[56 + i] + par[10]*par[3]*(1 + par[3])*par[16]*sin(x[4])*x[64 + i])/(1 + par[3]);
            f[40 + i] = x[48 + i];
            f[48 + i] = (cos(x[4])*(par[7]*par[7]*x[8 + i] + 2*par[4]*x[16 + i] - (par[10]*(1 + par[3])*par[8]*par[8] + x[2] + 2*par[5]*x[3] - par[15]*x[6])*x[40 + i]) - 
                         sin(x[4])*(x[24 + i] + 2*par[5]*x[32 + i] + (par[7]*par[7]*x[0] + 2*par[4]*x[1])*x[40 + i] - par[15]*x[56 + i]) - 
                         par[10]*(1 + par[3])*(par[9]*par[9]*x[40 + i] + 2*par[6]*x[48 + i] - par[16]*x[64 + i]))/par[10];
            f[56 + i] = -(par[12]*x[32 + i]) - par[11]*x[56 + i];
            f[64 + i] = -(par[14]*x[48 + i]) - par[13]*x[64 + i];
        }
    }
    else {
        error();
    }
}

void pend_oscillator_wout_pend_EH(int dim, double *x, double t, double *par, double *f) {
    // OMEGA   = par[0]   |   zeta_z    = par[5]          |   x[0] = x       
    // gamma   = par[1]   |   OMEGA_s   = par[6]          |   x[1] = dx/dt   
    // mu      = par[2]   |   varphi    = par[7]          |   x[2] = z       
    // rho     = par[3]   |   kappa     = par[8]          |   x[3] = dz/dt   
    // zeta_x  = par[4]   |   chi       = par[9]          |   x[4] = v        
    const double pi = 4 * atan(1);
    
    if (dim == 5) {
        f[0] = x[1];
        f[1] = (1/(1 + par[3]))*(-2*par[4]*x[1] - par[6]*par[6]*x[0]) + par[1]*par[0]*par[0]*sin(par[0]*t)*sin((pi/180)*par[2]); 
        f[2] = x[3];
        f[3] = (1/(1 + par[3]))*(-2*par[5]*x[3] - x[2] + par[9]*x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*cos((pi/180)*par[2]);
        f[4] = -par[8]*x[3] - par[7]*x[4];
    }
    else if (dim == 30) {
        f[0] = x[1];
        f[1] = (1/(1 + par[3]))*(-2*par[4]*x[1] - par[6]*par[6]*x[0]) + par[1]*par[0]*par[0]*sin(par[0]*t)*sin((pi/180)*par[2]); 
        f[2] = x[3];
        f[3] = (1/(1 + par[3]))*(-2*par[5]*x[3] - x[2] + par[9]*x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*cos((pi/180)*par[2]);
        f[4] = -par[8]*x[3] - par[7]*x[4];
        for (int i = 0; i < 5; i ++) {
            f[5 + i] = x[10 + i];
            f[10 + i] = -(1/(1 + par[3]))*(par[6]*par[6]*x[5+i] + 2*par[4]*x[10+i]);
            f[15 + i] = x[20 + i];
            f[20 + i] = -(1/(1 + par[3]))*(x[15+i] + 2*par[5]*x[20+i] - par[9]*x[25+i]);
            f[25 + i] = -par[8]*x[20+i] - par[7]*x[25+i];
        }
    }
    else {
        error();
    }
}

void duffing_2DoF_EH_old(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA   = par[0]   |   alpha_1 = par[5]   |   chi_1    = par[10]   |   kappa_2 = par[15]   |   x1  = x[0] | v2 = x[5]
       gamma   = par[1]   |   alpha_2 = par[6]   |   chi_2    = par[11]   |                       |   dx1 = x[1] |    
       rho     = par[2]   |   beta_1  = par[7]   |   varphi_1 = par[12]   |                       |   x2  = x[2] |
       zeta_1  = par[3]   |   beta_2  = par[8]   |   varphi_2 = par[13]   |                       |   dx2 = x[3] |
       zeta_2  = par[4]   |   OMEGA_s = par[9]   |   kappa_1  = par[14]   |                       |   v1  = x[4] |           */
    if (dim == 6) { 
        f[0] = x[1];
        f[1] = par[1]*par[0]*par[0]*sin(par[0] * t) - 2*par[3]*x[1] + 2*par[4]*(x[3] - x[1]) - (1 + par[5])*x[0] - par[7]*x[0]*x[0]*x[0]
               + par[9]*par[9]*par[2]*(x[2] - x[0]) + par[10]*x[4];
        f[2] = x[3];
        f[3] = par[1]*par[0]*par[0]*sin(par[0] * t) - (1/par[2])*(2*par[4]*(x[3] - x[1]) + par[6]*x[2] + par[8]*x[2]*x[2]*x[2] - par[11]*x[5])
               - par[9]*par[9]*(x[2] - x[0]);
        f[4] = -par[12]*x[4] - par[14]*x[1];
        f[5] = -par[13]*x[5] - par[15]*(x[3] - x[1]); 
    }
    else if (dim == 42) {
        f[0] = x[1];
        f[1] = par[1]*par[0]*par[0]*sin(par[0] * t) - 2*par[3]*x[1] + 2*par[4]*(x[3] - x[1]) - (1 + par[5])*x[0] - par[7]*x[0]*x[0]*x[0]
               + par[9]*par[9]*par[2]*(x[2] - x[0]) + par[10]*x[4];
        f[2] = x[3];
        f[3] = par[1]*par[0]*par[0]*sin(par[0] * t) - (1/par[2])*(2*par[4]*(x[3] - x[1]) + par[6]*x[2] + par[8]*x[2]*x[2]*x[2] - par[11]*x[5])
               - par[9]*par[9]*(x[2] - x[0]);
        f[4] = -par[12]*x[4] - par[14]*x[1];
        f[5] = -par[13]*x[5] - par[15]*(x[3] - x[1]); 
        for (int i = 0; i < 6; i ++) {
            f[6 + i] = x[12 + i];
            f[12 + i] = -((1 + par[5] + par[9]*par[9]*par[2] + 3*par[7]*x[0]*x[0])*x[6 + i]) - 2*(par[3] + par[4])*x[12 + i] + 
                        par[9]*par[9]*par[2]*x[18 + i] + 2*par[4]*x[24 + i] + par[10]*x[30 + i];
            f[18 + i] = x[24 + i];
            f[24 + i] = (par[9]*par[9]*par[2]*x[6 + i] + 2*par[4]*x[12 + i] - 
                        (par[6] + par[9]*par[9]*par[2] + 3*par[8]*x[2]*x[2])*x[18 + i] - 2*par[4]*x[24 + i] + par[11]*x[36 + i]
                        )/par[2];
            f[30 + i] = -par[14]*x[12 + i] - par[12]*x[30 + i];
            f[36 + i] = -par[15]*(x[24 + i] - x[12 + i]) - par[13]*x[36 + i];
        }
    }
    else {
        error();
    }
}

void duffing_2DoF_EH(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA   = par[0]   |   alpha_1 = par[5]   |   chi_1    = par[10]   |   kappa_2 = par[15]   |   x1  = x[0] | v2 = x[5]
       gamma   = par[1]   |   alpha_2 = par[6]   |   chi_2    = par[11]   |                       |   dx1 = x[1] |    
       rho     = par[2]   |   beta_1  = par[7]   |   varphi_1 = par[12]   |                       |   x2  = x[2] |
       zeta_1  = par[3]   |   beta_2  = par[8]   |   varphi_2 = par[13]   |                       |   dx2 = x[3] |
       zeta_2  = par[4]   |   OMEGA_s = par[9]   |   kappa_1  = par[14]   |                       |   v1  = x[4] |           */
    double ddxb = par[1]*par[0]*par[0]*sin(par[0] * t);
    if (dim == 6) { 
        f[0] = x[1]; 
        f[1] = ddxb - 2*par[3]*x[1] + 2*par[4]*(x[3] - x[1]) - (1 + par[5])*x[0] - par[7]*x[0]*x[0]*x[0]
               + par[9]*par[9]*par[2]*(x[2] - x[0]) + par[10]*x[4] - par[11]*x[5];
        f[2] = x[3];
        f[3] = ddxb - (1/par[2])*(2*par[4]*(x[3] - x[1]) + par[6]*x[2] + par[8]*x[2]*x[2]*x[2] - par[11]*x[5])
               - par[9]*par[9]*(x[2] - x[0]);
        f[4] = -par[12]*x[4] - par[14]*x[1];
        f[5] = -par[13]*x[5] - par[15]*(x[3] - x[1]); 
    }
    else if (dim == 42) {
        f[0] = x[1]; 
        f[1] = ddxb - 2*par[3]*x[1] + 2*par[4]*(x[3] - x[1]) - (1 + par[5])*x[0] - par[7]*x[0]*x[0]*x[0]
               + par[9]*par[9]*par[2]*(x[2] - x[0]) + par[10]*x[4] - par[11]*x[5];
        f[2] = x[3];
        f[3] = ddxb - (1/par[2])*(2*par[4]*(x[3] - x[1]) + par[6]*x[2] + par[8]*x[2]*x[2]*x[2] - par[11]*x[5])
               - par[9]*par[9]*(x[2] - x[0]);
        f[4] = -par[12]*x[4] - par[14]*x[1];
        f[5] = -par[13]*x[5] - par[15]*(x[3] - x[1]); 
        for (int i = 0; i < 6; i ++) {
            f[6 + i] = x[12 + i];
            f[12 + i] = -(1 + par[5] + par[9]*par[9]*par[2] + 3*par[7]*x[0]*x[0])*x[6 + i] - 2*(par[3] + par[4])*x[12 + i] + 
                        par[9]*par[9]*par[2]*x[18 + i] + 2*par[4]*x[24 + i] + par[10]*x[30 + i] - par[11]*x[36 + i];
            f[18 + i] = x[24 + i];
            f[24 + i] = (1/par[2])*(par[9]*par[9]*par[2]*x[6 + i] + 2*par[4]*x[12 + i] - (par[6] + par[9]*par[9]*par[2] + 3*par[8]*x[2]*x[2])*x[18 + i] 
                                    - 2*par[4]*x[24 + i] + par[11]*x[36 + i]);
            f[30 + i] = -par[14]*x[12 + i] - par[12]*x[30 + i];
            f[36 + i] = -par[15]*(x[24 + i] - x[12 + i]) - par[13]*x[36 + i];
        }
    }
    else {
        error();
    }
}

void linear_2DoF_EH(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA   = par[0]   |   OMEGA_s  = par[5]   |   kappa_1 = par[10]   |   x1  = x[0] |  v2 = x[5]
       gamma   = par[1]   |   chi_1    = par[6]   |   kappa_2 = par[11]   |   dx1 = x[1] |    
       rho     = par[2]   |   chi_2    = par[7]   |                       |   x2  = x[2] |
       zeta_1  = par[3]   |   varphi_1 = par[8]   |                       |   dx2 = x[3] |
       zeta_2  = par[4]   |   varphi_2 = par[9]   |                       |   v1  = x[4] |         */
    double ddxb = par[1]*par[0]*par[0]*sin(par[0] * t);
    if (dim == 6) { 
        f[0] = x[1];
        f[1] = ddxb - 2*par[3]*x[1] + 2*par[4]*(x[3] - x[1]) - x[0] + par[2]*par[5]*par[5]*(x[2] - x[0]) + par[6]*x[4] - par[7]*x[5];
        f[2] = x[3];
        f[3] = ddxb - (2/par[2])*par[4]*(x[3] - x[1]) - par[5]*par[5]*(x[2] - x[0])+ par[7]*x[5];
        f[4] = -par[8]*x[4] - par[10]*x[1];
        f[5] = -par[9]*x[5] - par[11]*(x[3] - x[1]); 
    }
    else if (dim == 42) {
        f[0] = x[1];
        f[1] = ddxb - 2*par[3]*x[1] + 2*par[4]*(x[3] - x[1]) - x[0] + par[2]*par[5]*par[5]*(x[2] - x[0]) + par[6]*x[4] - par[7]*x[5];
        f[2] = x[3];
        f[3] = ddxb - (2/par[2])*par[4]*(x[3] - x[1]) - par[5]*par[5]*(x[2] - x[0])+ par[7]*x[5];
        f[4] = -par[8]*x[4] - par[10]*x[1];
        f[5] = -par[9]*x[5] - par[11]*(x[3] - x[1]); 
        for (int i = 0; i < 6; i ++) {
            f[6 + i] = x[12 + i];
            f[12 + i] = -(1 + par[5]*par[5]*par[2])*x[6 + i] - 2*(par[3] + par[4])*x[12 + i] + par[5]*par[5]*par[2]*x[18 + i] 
                        + 2*par[4]*x[24 + i] + par[6]*x[30 + i] - par[7]*x[36 + i];
            f[18 + i] = x[24 + i];
            f[24 + i] = (1/par[2])*(par[5]*par[5]*par[2]*x[6 + i] - 2*par[4]*(x[24 + i] - x[12 + i]) - par[5]*par[5]*par[2]*x[18 + i]
                                    + par[7]*x[36 + i]);
            f[30 + i] = -par[10]*x[12 + i] - par[8]*x[30 + i];
            f[36 + i] = -par[11]*(x[24 + i] - x[12 + i]) - par[9]*x[36 + i];
        }
    }
    else {
        error();
    }
}

/* Adeodato (2020) SMA Model*/

static double sma_strain(double displ, double t, double t0, double strain_0, double L0) {
    // If initial time, then get initial condition properties
    if (t == t0) {
         return strain_0;
    } 
    // Else, update properties based on system behavior
    else {
        return (displ/L0);
    }
}

static double sma_stress(double t, double t0, double E, double sigma_0, double strain) {
    // If initial time, then get initial condition properties
    if (t == t0) {
         return sigma_0;
    } 
    // Else, update properties based on system behavior
    else {
        return E*strain;
    }
}

static double loop_direction(double t, double t0, double dstrain) {
    // If initial time, then get initial condition properties
    if (t == t0) {
         return 1.0;
    } 
    // Else, update properties based on system behavior
    else {
        return (dstrain/fabs(dstrain));
    }
}

static double load_direction(double t, double t0, double strain) {
    // If initial time, then get initial condition properties
    if (t == t0) {
         return 1.0;
    } 
    // Else, update properties based on system behavior
    else {
        return (strain/fabs(strain));
    }
}

static void IC_for_parameters(double *par_0, double *par, double *par_ant) {
    (*par_ant) = (*par_0);
    (*par) = (*par_0);
}

static double update_parameter(double par) {
    return par;
}

void adeodato_sma_oscillator(int dim, double *x, double t, double *par, double *f) {
    /* System Parameters  |                                                              Shape Memory Properties                                                               |   
       -----------------  -   ----------------------------------------------------------------------------------------------------------------------------------------------   -      
       OMEGA   = par[0]   |   sigma_0    = par[6]   |   beta_0    = par[12]   |   alpha  = par[18]  |  s_2      = par[24]   |   load_ant = par[30]   |   T         = par[36]   |   
       gamma   = par[1]   |   sigma_ant  = par[7]   |   beta_ant  = par[13]   |   Ms     = par[19]  |  Ca       = par[25]   |   load     = par[31]   |   t_ant     = par[37]   |      
       c       = par[2]   |   sigma      = par[8]   |   beta      = par[14]   |   Mf     = par[20]  |  Cm       = par[26]   |   Area     = par[32]   |   E_ant     = par[38]   |           
       k       = par[3]   |   strain_0   = par[9]   |   E_0       = par[15]   |   As     = par[21]  |  strain_r = par[27]   |   L_0      = par[33]   |   alpha_ant = par[39]   |
       m       = par[4]   |   strain_ant = par[10]  |   E         = par[16]   |   Af     = par[22]  |  dir_ant  = par[28]   |   Ea       = par[34]   |                         |
       t_0     = par[5]   |   strain     = par[11]  |   alpha_0   = par[17]   |   s_1    = par[23]  |  dir      = par[29]   |   Em       = par[35]   |                         |   */
    
    // Make parameter_ant = parameter_0 if t = t_0
    if (t == par[5]) {
        par[37] = update_parameter(par[5]);              // t_ant = t_0
        IC_for_parameters(&par[6], &par[8], &par[7]);    // Initial conditions for sigma
        IC_for_parameters(&par[9], &par[11], &par[10]);  // Initial conditions for strain
        IC_for_parameters(&par[12], &par[14], &par[13]); // Initial conditions for beta
        IC_for_parameters(&par[15], &par[16], &par[38]); // Initial conditions for E
        IC_for_parameters(&par[17], &par[18], &par[39]); // Initial conditions for alpha
        par[28] = update_parameter(par[29]);  // dir_ant = dir
        par[31] = update_parameter(par[30]);  // load_and = load
    }
    // if t bigger than t_ant, update properties
    if (t > par[37]) {
        // Update ant parameters
        par[7] = update_parameter(par[8]);    // sigma_ant = sigma
        par[10] = update_parameter(par[11]);  // strain_ant = strain
        par[13] = update_parameter(par[14]);  // beta_ant = beta
        par[38] = update_parameter(par[16]);  // E_ant = E
        par[39] = update_parameter(par[18]);  // alpha_ant = alpha
        par[28] = update_parameter(par[29]);  // dir_ant = dir
        par[30] = update_parameter(par[31]);  // load_ant = load
        // Compute current sma strain
        par[11] = sma_strain(x[0], t, par[5], par[9], par[33]);    // Strain
        // Compute strain derivative
        double dstrain = par[11] - par[10];
        // Compute directions of steps in the loop and direction of the load
        par[29] = loop_direction(t, par[5], dstrain);    // dir
        par[31] = load_direction(t, par[5], par[11]);    // load
        // Update parameters if load or dir changes signal
        if ((par[29] != par[28]) || (par[31] != par[30])) {
            printf("t = %lf | dir_ant = %lf | dir = %lf\n", t, par[28], par[29]);
            printf("t = %lf | load_ant = %lf | load = %lf\n", t, par[30], par[31]);
            par[6] = update_parameter(par[7]);      // sigma_0 = sigma_ant
            par[9] = update_parameter(par[10]);     // strain_0 = strain_ant
            par[12] = update_parameter(par[13]);    // beta_0 = beta_ant
            par[15] = update_parameter(par[38]);    // E_0 = E_ant
            par[17] = update_parameter(par[39]);    // alpha_0 = alpha_ant
        }
        // Compute current sma stress
        par[8] = sma_stress(t, par[5], par[16], par[6], par[11]);  // Sigma
        
    }
    // System of Equations
    if (dim == 2) {
        f[0] = x[1];
        f[1] = (1/par[4])*(par[1]*sin(par[0] * t) - par[2]*x[1] - par[3]*x[0] - par[8]*par[32]);    
    } 
    else if (dim == 4) {
        f[0] = x[1];
        f[1] = (1/par[4])*(par[1]*sin(par[0] * t) - par[2]*x[1] - par[3]*x[0] - par[8]*par[32]);
        for (int i = 0; i < 2; i++) {
            f[2 + i] = 0;
            f[4 + i] = 0;
        }
    }
    else {
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
    }
}

// Not Implemented
void duffing_cldyn(int dim, double *x, double t, double *par, double *f) {
    // OMEGA = par[0]
    // gamma = par[1]
    // zeta = par[2]
    // alpha = par[3]
    // beta = par[4]
    
    if (dim == 2) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0];    
    } 
    else if (dim == 6) {
        f[0] = x[1];
        f[1] = par[1]*sin(par[0] * t) - 2*par[2]*x[1] - par[3]*x[0] - par[4]*x[0]*x[0]*x[0];
        for (int i = 0; i < 2; i++) {
            f[2 + i] = x[4 + i];
            f[4 + i] = par[1]*sin(par[0] * t) - 2*par[2]*x[4 + i] - par[3]*x[2 + i] - par[4]*x[2 + i]*x[2 + i]*x[2 + i];
        }
    }
    else {
        error();
    }    
}
