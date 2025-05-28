#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif

// #define setm->omega_aso 5.0E-3
// #define setm->mass_aso 20000.0

// int type_aso;
// double setm->A_aso[9], setm->beta_aso[3], setm->D_aso[3], setm->C_aso[3], setm->R0_aso[3];
// double setm->R00_aso[9], setm->alpha_aso[9];
// double setm->Re_aso;

// void readinp_aso(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &type_aso);
//     *Ndof1 = 1;
//     *Ndof2 = 1;
//     *Nstate = 3;
// }

void parameter_aso(double *mass, struct set_host *setm) {
    int i, j;

    setm->mass_aso = (74.922*15.999)/(74.922+15.999) * amu_2_au;
    setm->omega_aso = 1160/au_2_wn;

    mass[0] = setm->mass_aso;

    // double S0[0] = { 1.01597749,  2.96496059, -0.00757351,0.34817623};
    // double S1[0] = { 1.01046726 , 3.27798996, 0.15125878, 0.15498669};
    // double T1[0] = { 1.05112708 , 3.2378816,  0.13893564, 0.14709944};

    setm->Re_aso = 2.5;

    setm->A_aso[0] = 1.01597749; // S0
    setm->A_aso[1] = 1.01046726; // S1
    setm->A_aso[2] = 1.05112708; // T1

    setm->R0_aso[0] = 2.96496059; // S0
    setm->R0_aso[1] = 3.27798996; // S1
    setm->R0_aso[2] = 3.2378816; // T1

    setm->S_aso[0] = -0.00757351; // S0
    setm->S_aso[1] = 0.15125878; // S1
    setm->S_aso[2] = 0.13893564; // T1

    setm->d_aso[0] = 0.34817623; // S0
    setm->d_aso[1] = 0.15498669; // S1
    setm->d_aso[2] = 0.14709944; // T1 

    setm->c_aso = 2.5;
    setm->b_aso = 10; 
    setm->rnac_aso = 4.9;




    setm->soc_aso = (double complex *)malloc(5 * 5 * sizeof(double complex));
    memset(setm->soc_aso, 0, 5 * 5 * sizeof(double complex)); 

    setm->soc_aso[0 * 5 + 2] = (2500.0 + 0.0 * I)/au_2_wn ; // S0 - T1, ms=-1
    setm->soc_aso[0 * 5 + 3] = (0.0 + 2000.0 * I)/au_2_wn ; // S0 - T1, ms=0
    setm->soc_aso[0 * 5 + 4] = (2500.0 + 0.0 * I)/au_2_wn ; // S0 - T1, ms=1
    setm->soc_aso[1 * 5 + 2] = (0.0 - 1500.0 * I)/au_2_wn ; // S1 - T1, ms=-1
    setm->soc_aso[1 * 5 + 3] = (0.0 + 0.0 * I)/au_2_wn ; // S1 - T1, ms=0
    setm->soc_aso[1 * 5 + 4] = (0.0 + 1500.0 * I)/au_2_wn ; // S1 - T1, ms=1
    

    setm->soc_aso[ 2 * 5 + 0 ] = (2500.0 + 0.0 * I)/au_2_wn ; 
    setm->soc_aso[ 3 * 5 + 0 ] = (0.0 - 2000.0 * I)/au_2_wn ; 
    setm->soc_aso[ 4 * 5 + 0 ] = (2500.0 + 0.0 * I)/au_2_wn ; 
    setm->soc_aso[ 2 * 5 + 1 ] = (0.0 + 1500.0 * I)/au_2_wn ; 
    setm->soc_aso[ 3 * 5 + 1 ] = (0.0 + 0.0 * I)/au_2_wn ; 
    setm->soc_aso[ 4 * 5 + 1 ] = (0.0 - 1500.0 * I)/au_2_wn ; 
    


}

void sample_aso(double *P, double *R, struct set_host *setm) {
    double x2, gamma_aso;
    gamma_aso = setm->mass_aso * setm->omega_aso;

    // Assuming box_muller is defined elsewhere
    box_muller(&P[0], &x2, sqrt(gamma_aso / 2), 0.0);
    box_muller(&R[0], &x2, sqrt(1 / (gamma_aso * 2)), setm->Re_aso);
}

void V_aso(double *R, double complex *H, struct set_host *setm) {
    int i, j;

    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            if (i == j) {
                if (i < 2) {
                    // Diagonal elements for S0, S1
                    H[i * 5 + i] = setm->d_aso[i] * pow(1 - exp(-setm->A_aso[i] * (R[0] - setm->R0_aso[i])), 2) + setm->S_aso[i];
                } else {
                    // Diagonal elements for T1
                    H[i * 5 + i] = setm->d_aso[2] * pow(1 - exp(-setm->A_aso[2] * (R[0] - setm->R0_aso[2])), 2) + setm->S_aso[2]; 
                }
                
            } else {
                H[i * 5 + j] = setm->soc_aso[i * 5 + j];
            }
        }
    }

}

void dV_aso(double *R, double complex *dH, struct set_host *setm) {
    int i, j;

    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            if (i == j) {
                if (i < 2) {
                    // Diagonal elements for S0, S1
                    dH[i * 5 + i] = setm->d_aso[i] * 2 * (1 - exp(-setm->A_aso[i] * (R[0] - setm->R0_aso[i]))) * setm->A_aso[i] * exp(-setm->A_aso[i] * (R[0] - setm->R0_aso[i]));
                } else {
                    // Diagonal elements for T1
                    dH[i * 5 + i] = setm->d_aso[2] * 2 * (1 - exp(-setm->A_aso[2] * (R[0] - setm->R0_aso[2]))) * setm->A_aso[2] * exp(-setm->A_aso[2] * (R[0] - setm->R0_aso[2]));
                }
               
            } else {
                dH[i * 5 + j] = 0.0 + I * 0.0; 
            }
        }
    }

}



void nac_aso(double *R, double complex *nac, struct set_host *setm) {

    memset(nac, 0, 5 * 5 * sizeof(double)); // Initialize nac to zero
    nac[0 * 5 + 1] = setm->c_aso * exp(-setm->b_aso * (R[0] - setm->rnac_aso) * (R[0] - setm->rnac_aso)); // S0 - S1
    nac[1 * 5 + 0] = -nac[0 * 5 + 1]; // S1 - S0
    
}
