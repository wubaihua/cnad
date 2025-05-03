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
#include "lvcm.h"

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


void parameter_retinal(double *mass, struct set_host *setm) {
    
    parameter_LVCM(mass, setm);

    // mode 0-22: bath mode
    // mode 23: qc, mode 24: theta

    setm->H_ele_lvcm[3] = 2.48 / au_2_eV ;
    setm->W1_retinal = 3.6 / au_2_eV ;
    setm->W2_retinal = 1.09 / au_2_eV ;

    mass[24] = 56198.347;


    double w[25] = {792.8,842.8,866.2,882.4,970.3,976.0,997.0,1017.1,1089.6,1189.0,
                    1214.7,1238.1,1267.9,1317.0,1359.0,1389.0,1428.4,1434.9,1451.8,
                    1572.8,1612.1,1629.2,1659.1, 1532.0, 0.0};
            
    
    double c[25] = {0.175,0.2,0.175,0.225,0.55,0.3,0.33,0.45,0.125,0.175,
                    0.44,0.5,0.475,0.238,0.25,0.25,0.25,0.225,0.225,0.25,
                    0.225,0.125,0.225,0.0,0.0};

    for (int i = 0; i < 25; i++){
        setm->omega_lvcm[i] = w[i] / au_2_wn;
        setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + i] = c[i] * w[i]/ au_2_wn;
    }
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.1 / au_2_eV;

    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.19 / au_2_eV;

    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.19 / au_2_eV;

    
     
}

void sample_retinal(double *P, double *R, struct set_host *setm) {
    int j;
    double x2,sg;

    sample_LVCM(P, R, setm);
    sg=0.15228;
    box_muller(&P[24], &x2, 1.0 / sg , 0.0);
    // P[24] = 0.0;
    box_muller(&R[24], &x2, sg / 2.0, 0.0);

    

}

void V_retinal(double *R, double complex *H, int forcetype, struct set_host *setm) {
    
    V_LVCM(R, H, forcetype, setm);

    H[0] += setm->W1_retinal/2 * (1.0 - cos(R[24]));
    H[3] += -setm->W2_retinal/2 * (1.0 - cos(R[24]));

    
    
}

void dV_retinal(double *R, double complex *dH, int forcetype, struct set_host *setm) {
    
    dV_LVCM(R, dH, forcetype, setm);

    dH[0 * setm->Nstate_lvcm * setm->N_mode_lvcm +0 * setm->N_mode_lvcm+ 24] += setm->W1_retinal/2 * sin(R[24]);
    dH[1 * setm->Nstate_lvcm * setm->N_mode_lvcm +1 * setm->N_mode_lvcm+ 24] += -setm->W2_retinal/2 * sin(R[24]);

    

}

void nucforce_retinal(double *R, double *nf, struct set_host *setm) {

    nucforce_LVCM(R, nf, setm);

}



