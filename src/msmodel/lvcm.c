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

// int setm->N_mode_lvcm, setm->Nstate_lvcm;
// double setm->L_lvcm;
// double *setm->eps_lvcm, *setm->miu_lvcm, *setm->omega_lvcm, *setm->lambda_lvcm;

// void readinp_lvcm(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &setm->Nstate_lvcm);
//     fscanf(idinp, "%d", &setm->N_mode_lvcm);

//     *Ndof1 = 1;
//     *Ndof2 = setm->N_mode_lvcm;
//     *Nstate = setm->Nstate_lvcm;
// }


void parameter_LVCM(double *mass, struct set_host *setm) {
    
    int j;

    for (j = 0; j < setm->N_mode_lvcm; j++) {
        mass[j] = 1.0;
    }

    setm->H_ele_lvcm = (double *)malloc(setm->Nstate_lvcm * setm->Nstate_lvcm * sizeof(double));
    memset(setm->H_ele_lvcm, 0, setm->Nstate_lvcm * setm->Nstate_lvcm * sizeof(double));
    setm->c_lvcm = (double *)malloc(setm->Nstate_lvcm * setm->Nstate_lvcm * setm->N_mode_lvcm * sizeof(double));
    memset(setm->c_lvcm, 0, setm->Nstate_lvcm * setm->Nstate_lvcm * setm->N_mode_lvcm * sizeof(double));
    setm->omega_lvcm = (double *)malloc(setm->N_mode_lvcm * sizeof(double));
    memset(setm->omega_lvcm, 0, setm->N_mode_lvcm * sizeof(double));
    setm->R0_lvcm = (double *)malloc(setm->N_mode_lvcm * sizeof(double));
    memset(setm->R0_lvcm, 0, setm->N_mode_lvcm * sizeof(double));
    setm->P0_lvcm = (double *)malloc(setm->N_mode_lvcm * sizeof(double));
    memset(setm->P0_lvcm, 0, setm->N_mode_lvcm * sizeof(double));
    setm->alpha_lvcm = (double *)malloc(setm->N_mode_lvcm * sizeof(double));
    for (j = 0; j < setm->N_mode_lvcm; j++) {
        setm->alpha_lvcm[j] = 1.0;
    }
}

void sample_LVCM(double *P, double *R, struct set_host *setm) {
    int j;
    double x2;

    for (j = 0; j < setm->N_mode_lvcm; j++) {
        // Assuming box_muller is a function that generates Gaussian random numbers
        box_muller(&P[j], &x2, sqrt(0.5 * setm->omega_lvcm[j])/setm->alpha_lvcm[j], setm->P0_lvcm[j]);
        box_muller(&R[j], &x2, setm->alpha_lvcm[j] * sqrt(0.5 / setm->omega_lvcm[j]), setm->R0_lvcm[j]);
    }

    

}

void V_LVCM(double *R, double *H, int forcetype, struct set_host *setm) {
    int i, j, k;

    for (i = 0; i < setm->Nstate_lvcm * setm->Nstate_lvcm; i++) {
        H[i] = setm->H_ele_lvcm[i];
    }

    for (i = 0; i < setm->Nstate_lvcm; i++){
        for (j = 0; j < setm->Nstate_lvcm; j++){
            for (k = 0; k < setm->N_mode_lvcm; k++){
                H[i * setm->Nstate_lvcm + j] += setm->c_lvcm[i * setm->Nstate_lvcm * setm->N_mode_lvcm 
                                                             + j * setm->N_mode_lvcm + k] * sqrt(setm->omega_lvcm[k]) * R[k];
            }
        }
    }

    if (forcetype == 0) {
        double sum = 0;
        for (i = 0; i < setm->N_mode_lvcm; i++){
            sum += 0.5 * setm->omega_lvcm[i] * setm->omega_lvcm[i] * R[i] * R[i];
        }
        for (i = 0; i < setm->Nstate_lvcm; i++) {
            H[i * setm->Nstate_lvcm + i] += sum;
        }
    }

    
}

void dV_LVCM(double *R, double *dH, int forcetype, struct set_host *setm) {
    int i, j, k;

    for (i = 0; i < setm->Nstate_lvcm; i++){
        for (j = 0; j < setm->Nstate_lvcm; j++){
            for (k = 0; k < setm->N_mode_lvcm; k++){
                dH[i * setm->Nstate_lvcm * setm->N_mode_lvcm + j * setm->N_mode_lvcm + k] 
                = setm->c_lvcm[i * setm->Nstate_lvcm * setm->N_mode_lvcm + j * setm->N_mode_lvcm + k] 
                  * sqrt(setm->omega_lvcm[k]);
            }
        }
    }

   

    if (forcetype == 0) {
        for (i = 0; i < setm->Nstate_lvcm; i++) {
            for (j = 0; j < setm->N_mode_lvcm; j++) {
               
                dH[i * setm->Nstate_lvcm * setm->N_mode_lvcm +i * setm->N_mode_lvcm+ j] += setm->omega_lvcm[j] * setm->omega_lvcm[j] * R[j];
               
            }
        }   
    }
}

void nucforce_LVCM(double *R, double *nf, struct set_host *setm) {
    int j;

    for (j = 0; j < setm->N_mode_lvcm; j++) {
        nf[j] = setm->omega_lvcm[j] * setm->omega_lvcm[j] * R[j];
    }
}



