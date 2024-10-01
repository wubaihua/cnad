#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"
#include <string.h>
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif

// #define PI 3.141592653589793
// #define HBAR 1.0545718e-34

// int setm->N_mode_aic, setm->Nstate_aic;
// double setm->L_aic;
// double *setm->eps_aic, *setm->miu_aic, *setm->omega_aic, *setm->lambda_aic;

// void readinp_AIC(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &setm->Nstate_aic);
//     fscanf(idinp, "%d", &setm->N_mode_aic);

//     *Ndof1 = 1;
//     *Ndof2 = setm->N_mode_aic;
//     *Nstate = setm->Nstate_aic;
// }

void parameter_AIC(double *mass, struct set_host *setm) {
    double c, eps0;
    int j;

    for (j = 0; j < setm->N_mode_aic; j++) {
        mass[j] = 1.0;
    }

    setm->eps_aic = (double *)malloc(setm->Nstate_aic * sizeof(double));
    setm->miu_aic = (double *)malloc(setm->Nstate_aic * setm->Nstate_aic * sizeof(double));
    setm->omega_aic = (double *)malloc(setm->N_mode_aic * sizeof(double));
    setm->lambda_aic = (double *)malloc(setm->N_mode_aic * sizeof(double));

    eps0 = 1.0 / (4 * M_PI);
    c = 137.035999139;
    setm->L_aic = 2.362e5;

    for (j = 0; j < setm->N_mode_aic; j++) {
        setm->omega_aic[j] = (j + 1) * M_PI * c / setm->L_aic;
        setm->lambda_aic[j] = sqrt(2.0 / (eps0 * setm->L_aic)) * sin(0.5 * setm->L_aic * setm->omega_aic[j] / c);
    }

    for (j = 0; j < setm->Nstate_aic * setm->Nstate_aic; j++) {
        setm->miu_aic[j] = 0.0;
    }

    if (setm->Nstate_aic == 2) {
        setm->eps_aic[0] = -0.6738;
        setm->eps_aic[1] = -0.2798;
        setm->miu_aic[1 * setm->Nstate_aic + 0] = 1.034;
        setm->miu_aic[0 * setm->Nstate_aic + 1] = 1.034;
    } else if (setm->Nstate_aic == 3) {
        setm->eps_aic[0] = -0.6738;
        setm->eps_aic[1] = -0.2798;
        setm->eps_aic[2] = -0.1547;
        setm->miu_aic[1 * setm->Nstate_aic + 0] = 1.034;
        setm->miu_aic[2 * setm->Nstate_aic + 1] = -2.536;
        setm->miu_aic[0 * setm->Nstate_aic + 1] = 1.034;
        setm->miu_aic[1 * setm->Nstate_aic + 2] = -2.536;
    }
}

void sample_AIC(double *P, double *R, struct set_host *setm) {
    int j;
    double x2;

    for (j = 0; j < setm->N_mode_aic; j++) {
        // Assuming box_muller is a function that generates Gaussian random numbers
        box_muller(&P[j], &x2, sqrt(0.5 * hbar * setm->omega_aic[j]), 0.0);
        box_muller(&R[j], &x2, sqrt(0.5 * hbar / setm->omega_aic[j]), 0.0);
    }

    

}

void V_AIC(double *R, double *H, int forcetype, struct set_host *setm) {
    int i, j;

    for (i = 0; i < setm->Nstate_aic * setm->Nstate_aic; i++) {
        H[i] = 0.0;
    }

    double sum = 0;
    for (i = 0; i < setm->N_mode_aic; i++){
        sum += setm->omega_aic[i] * setm->lambda_aic[i] * R[i];
    }

    for (i = 0; i < setm->Nstate_aic; i++) {
        for (j = 0; j < setm->Nstate_aic; j++) {
            if (i == j) {
                H[i * setm->Nstate_aic + i] = setm->eps_aic[i];
            } else {
                H[i * setm->Nstate_aic + j] = setm->miu_aic[i * setm->Nstate_aic + j] * sum;
            }
        }
    }

    if (forcetype == 0) {
        sum = 0;
        for (i = 0; i < setm->N_mode_aic; i++){
            sum += 0.5 * setm->omega_aic[i] * setm->omega_aic[i] * R[i] * R[i];
        }
        for (i = 0; i < setm->Nstate_aic; i++) {
            H[i * setm->Nstate_aic + i] += sum;
        }
    }

    
}

void dV_AIC(double *R, double *dH, int forcetype, struct set_host *setm) {
    int i, j, k;

    for (i = 0; i < setm->Nstate_aic * setm->Nstate_aic * setm->N_mode_aic; i++) {
        dH[i] = 0.0;
    }

    for (i = 0; i < setm->Nstate_aic; i++) {
        for (j = 0; j < setm->Nstate_aic; j++) {
            if (i != j) {
                for (k = 0; k < setm->N_mode_aic; k++) {
                    dH[(i * setm->Nstate_aic + j) * setm->N_mode_aic + k] = setm->miu_aic[i * setm->Nstate_aic + j] * setm->omega_aic[k] * setm->lambda_aic[k];
                }
            }
        }
    }

    if (forcetype == 0) {
        for (i = 0; i < setm->Nstate_aic; i++) {
            for (j = 0; j < setm->N_mode_aic; j++) {
               
                dH[i * setm->Nstate_aic * setm->N_mode_aic +i * setm->N_mode_aic+ j] += setm->omega_aic[j] * setm->omega_aic[j] * R[j];
               
            }
        }   
    }
}

void nucforce_AIC(double *R, double *nf, struct set_host *setm) {
    int j;

    for (j = 0; j < setm->N_mode_aic; j++) {
        nf[j] = setm->omega_aic[j] * setm->omega_aic[j] * R[j];
    }
}



