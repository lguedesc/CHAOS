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
// 2 Nonlinear DoF EH
// Oscillator-Pendulum Energy Harvester

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
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
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
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
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
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
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
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
    }    
}

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
            f[3 + i] = 0;       // Add Linearized Eqs Later
            f[6 + i] = 0;       // Add Linearized Eqs Later
            f[9 + i] = 0;       // Add Linearized Eqs Later
        }
    }
    else {
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
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
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
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
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
    }    
}

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
        printf("Wrong dimension (dim) or (ndim) allocated for system of equations\n");
        exit(1);
    }    
}

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
}

void pend_oscillator_EH(int dim, double *x, double t, double *par, double *f) {
    /* OMEGA   = par[0]   |   zeta_z    = par[5]   |   chi    = par[10]               |   x[0] = x       |   x[5] = dphi/dt
       gamma   = par[1]   |   zeta_t    = par[6]   |   varphi = par[11]               |   x[1] = dx/dt   |   x[6] = v
       mu      = par[2]   |   OMEGA_s   = par[7]   |   kappa  = par[12]               |   x[2] = z       |
       rho     = par[3]   |   OMEGA_phi = par[8]   |                                  |   x[3] = dz/dt   |
       zeta_x  = par[4]   |   OMEGA_t   = par[9]   |                                  |   x[4] = phi     |                   */       
    if (dim == 7) {
        f[0] = x[1];
        f[1] = (1/(1 + par[3]))*(-(1 + par[3]*cos(x[4])*cos(x[4]))*(2*par[4]*x[1] + par[7]*par[7]*x[0]) + (par[3]/2)*(2*par[5]*x[3] + x[2] - par[10]*x[6])*sin(2*x[4]) + par[3]*x[5]*x[5]*sin(x[4]))
                + (2*par[6]*x[5] + par[9]*par[9]*x[4])*par[3]*cos(x[4]) + (par[3]/2)*par[8]*par[8]*sin(2*x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*sin(par[2]); 
        f[2] = x[3];
        f[3] = (1/(1 + par[3]))*((1 + par[3]*sin(x[4])*sin(x[4]))*(-2*par[5]*x[3] - x[2] + par[10]*x[6]) + (par[3]/2)*(2*par[4]*x[1] + par[7]*par[7]*x[0])*sin(2*x[4]) + par[3]*x[5]*x[5]*cos(x[4]))
                - (2*par[6]*x[5] + par[9]*par[9]*x[4])*par[3]*sin(x[4]) - par[3]*par[8]*par[8]*sin(x[4])*sin(x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*cos(par[2]);
        f[4] = x[5];
        f[5] = -(1 + par[3])*(2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4])) + (2*par[4]*x[1] + par[7]*par[7]*x[0])*cos(x[4]) - (2*par[5]*x[3] + x[2] - par[10]*x[6])*sin(x[4]);
        f[6] = -par[12]*x[3] - par[11]*x[6];
    }
    else if (dim == 56) {
        f[0] = x[1];
        f[1] = (1/(1 + par[3]))*(-(1 + par[3]*cos(x[4])*cos(x[4]))*(2*par[4]*x[1] + par[7]*par[7]*x[0]) + (par[3]/2)*(2*par[5]*x[3] + x[2] - par[10]*x[6])*sin(2*x[4]) + par[3]*x[5]*x[5]*sin(x[4]))
                + (2*par[6]*x[5] + par[9]*par[9]*x[4])*par[3]*cos(x[4]) + (par[3]/2)*par[8]*par[8]*sin(2*x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*sin(par[2]); 
        f[2] = x[3];
        f[3] = (1/(1 + par[3]))*((1 + par[3]*sin(x[4])*sin(x[4]))*(-2*par[5]*x[3] - x[2] + par[10]*x[6]) + (par[3]/2)*(2*par[4]*x[1] + par[7]*par[7]*x[0])*sin(2*x[4]) + par[3]*x[5]*x[5]*cos(x[4]))
                - (2*par[6]*x[5] + par[9]*par[9]*x[4])*par[3]*sin(x[4]) - par[3]*par[8]*par[8]*sin(x[4])*sin(x[4]) + par[1]*par[0]*par[0]*sin(par[0]*t)*cos(par[2]);
        f[4] = x[5];
        f[5] = -(1 + par[3])*(2*par[6]*x[5] + par[9]*par[9]*x[4] + par[8]*par[8]*sin(x[4])) + (2*par[4]*x[1] + par[7]*par[7]*x[0])*cos(x[4]) - (2*par[5]*x[3] + x[2] - par[10]*x[6])*sin(x[4]);
        f[6] = -par[12]*x[3] - par[11]*x[6];
        for (int i = 0; i < 7; i ++) {
            f[7 + i] = x[14 + i];
            f[14 + i] = 0;
            f[21 + i] = x[28 + i];
            f[28 + i] = 0;
            f[35 + i] = x[42 + i];
            f[42 + i] = 0;
            f[49 + i] = 0;
        }
    }
}