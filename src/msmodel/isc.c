#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif

// #define setm->omega_isc 5.0E-3
// #define setm->mass_isc 20000.0

// int type_isc;
// double setm->A_isc[9], setm->beta_isc[3], setm->D_isc[3], setm->C_isc[3], setm->R0_isc[3];
// double setm->R00_isc[9], setm->alpha_isc[9];
// double setm->Re_isc;

// void readinp_isc(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &type_isc);
//     *Ndof1 = 1;
//     *Ndof2 = 1;
//     *Nstate = 3;
// }

void parameter_isc(double *mass, struct set_host *setm) {
    int i, j;

   
    if (setm->type_isc == 1) {
        mass[0] = 20000.0; // in a.u.
        setm->a1_isc = 0.03452;
        setm->a2_isc = 0.5;
        setm->alpha1_isc = 0.35;
        setm->alpha2_isc = 0.25;
        setm->dE_isc = 0.04; 
        setm->c0_isc = 0.001 + 0.0 * I;
        setm->c1_isc = 0.0005 + 0.0005 * I; 
    }

    
 
    


}

void sample_isc(double *P, double *R, struct set_host *setm) {
    double x2, gamma_isc;


    if (setm->type_isc == 1) {
        gamma_isc = 20.0;
        box_muller(&P[0], &x2, sqrt(gamma_isc / 2), 20.0);
        box_muller(&R[0], &x2, sqrt(1 / (gamma_isc * 2)), 5.0);
    }
}

void V_isc(double *R, double complex *H, struct set_host *setm) {
    int i, j;


    if (setm->type_isc == 1) {

        memset(H, 0, 4 * 4 * sizeof(double complex)); // Initialize H to zero
        H[0 * 4 + 0] = setm->a1_isc * exp(-setm->alpha1_isc * R[0]) + setm->dE_isc; // S1
        H[1 * 4 + 1] = setm->a2_isc * exp(-setm->alpha2_isc * R[0]); // T1
        H[2 * 4 + 2] = H[1 * 4 + 1]; // T1
        H[3 * 4 + 3] = H[1 * 4 + 1]; // T1

        double sigma = 0.0;
        if (R[0] < setm->rs_isc - setm->drs_isc / 2) {
            sigma = 1.0;
        } else if (R[0] > setm->rs_isc + setm->drs_isc / 2) {
            sigma = -1.0;
        } else {
            sigma = 4.0 * pow((R[0] - setm->rs_isc) / setm->drs_isc, 3) - 3.0 * (R[0] - setm->rs_isc) / setm->drs_isc;
        }

        double complex z = setm->c1_isc * sigma;
        double complex b = setm->c0_isc * sigma;

        H[0 * 4 + 1] = z; 
        H[1 * 4 + 0] = conj(z);
        H[0 * 4 + 2] = I * b; 
        H[2 * 4 + 0] = -I * conj(b);
        H[0 * 4 + 3] = conj(z);
        H[3 * 4 + 0] = z;

    }

    

}

void dV_isc(double *R, double complex *dH, struct set_host *setm) {
    int i, j;

    if (setm->type_isc == 1) {

        memset(dH, 0, 4 * 4 * sizeof(double complex)); // Initialize H to zero
        dH[0 * 4 + 0] = -1.0 * setm->alpha1_isc * setm->a1_isc * exp(-setm->alpha1_isc * R[0]); 
        dH[1 * 4 + 1] = -1.0 * setm->alpha2_isc * setm->a2_isc * exp(-setm->alpha2_isc * R[0]); // T1
        dH[2 * 4 + 2] = dH[1 * 4 + 1]; // T1
        dH[3 * 4 + 3] = dH[1 * 4 + 1]; // T1

        double sigma = 0.0;
        if (R[0] < setm->rs_isc - setm->drs_isc / 2) {
            sigma = 0.0;
        } else if (R[0] > setm->rs_isc + setm->drs_isc / 2) {
            sigma = 0.0;
        } else {
            sigma = 12.0 * pow((R[0] - setm->rs_isc) / setm->drs_isc, 2) / setm->drs_isc - 3.0 / setm->drs_isc;
        }

        double complex z = setm->c1_isc * sigma;
        double complex b = setm->c0_isc * sigma;

        dH[0 * 4 + 1] = z; 
        dH[1 * 4 + 0] = conj(z);
        dH[0 * 4 + 2] = I * b; 
        dH[2 * 4 + 0] = -I * conj(b);
        dH[0 * 4 + 3] = conj(z);
        dH[3 * 4 + 0] = z;

    }

}



void nac_isc(double *R, double complex *nac, struct set_host *setm) {

    if (setm->type_isc == 1) {
        memset(nac, 0, 4 * 4 * sizeof(double)); // Initialize nac to zero
    }
    
    
}
