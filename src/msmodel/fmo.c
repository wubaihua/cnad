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
#include "sem.h"

// double lambda_SEM, omega_c_SEM;
// int N_bath_SEM;

// void readinp_FMO(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%lf", &lambda_SEM);
//     lambda_SEM = lambda_SEM * AU_2_PS * 10 / 53.09;
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%lf", &omega_c_SEM);
//     omega_c_SEM = omega_c_SEM * AU_2_PS * 10 / 53.09;
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%d", &N_bath_SEM);

//     *Nstate = NSTATE_SEM;
//     *Ndof1 = NSTATE_SEM;
//     *Ndof2 = N_bath_SEM;
// }

void parameter_FMO(double *mass, struct set_host *setm) {
    
   

    // 请根据实际情况实现parameter_SEM函数
    

    double H_ele[49] = {
        12410.0,     -87.7,      5.5,       -5.9,     6.7,       -13.7,        -9.9,
        -87.7,       12530.0,    30.8,      8.2,      0.7,        11.8,        4.3,
        5.5,         30.8,       12210.0,   -53.5,    -2.2,       -9.6,        6.0,
        -5.9,        8.2,        -53.5,     12320.0,  -70.7,      -17.0,       -63.3,
        6.7,         0.7,        -2.2,      -70.7,    12480.0,   81.1,       -1.3,
        -13.7,       11.8,       -9.6,      -17.0,    81.1,      12630.0,    39.7,
        -9.9,        4.3,        6.0,       -63.3,    -1.3,       39.7,      12440.0
    };


    setm->H_ele_SEM = (double *)malloc(setm->Nstate_SEM * setm->Nstate_SEM * sizeof(double));

    for (int i = 0; i < setm->Nstate_SEM * setm->Nstate_SEM; i++) {
        setm->H_ele_SEM[i] = H_ele[i] / au_2_wn;
    }


    for (int i = 0; i < setm->Nstate_SEM * setm->N_bath_SEM; i++) {
        mass[i] = 1.0;
    }

    setm->c_SEM = (double *)malloc(setm->N_bath_SEM * sizeof(double));
    setm->omega_SEM = (double *)malloc(setm->N_bath_SEM * sizeof(double));

    for (int j = 0; j < setm->N_bath_SEM; j++) {
        setm->omega_SEM[j] = setm->omega_c_SEM * tan(0.5 * M_PI * (1.0 - (double)(j + 1) / (setm->N_bath_SEM + 1)));
        setm->c_SEM[j] = setm->omega_SEM[j] * sqrt(setm->lambda_SEM * 2 / (setm->N_bath_SEM + 1));
    }



}

void sample_FMO(double *P, double *R, double beta, struct set_host *setm) {
    // 请根据实际情况实现sample_SEM函数
    sample_SEM(P, R, beta, setm);
}

void V_FMO(double *R, double *H, int forcetype, struct set_host *setm) {
    // 请根据实际情况实现V_SEM函数
    V_SEM(R, H, forcetype, setm);
}

void dV_FMO(double *R, double *dH, int forcetype, struct set_host *setm) {
    // 请根据实际情况实现dV_SEM函数
    dV_SEM(R, dH, forcetype, setm);
}

void nucforce_FMO(double *R, double *nf, struct set_host *setm) {
    // 请根据实际情况实现nucforce_SEM函数
    nucforce_SEM(R, nf, setm);
}

