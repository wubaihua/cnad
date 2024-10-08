#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"  // Assuming you have a math.h header for functions like box_muller
// #include "constant.h"
#include <string.h>
#include <stdio.h>
// #include "cJSON.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif
#include "semdp.h"

// void readinp_FMOdp(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &type_cal);
//     fscanf(idinp, "%lf", &lambda_SEMdp);
//     lambda_SEMdp = lambda_SEMdp / au_2_wn;
//     fscanf(idinp, "%lf", &omega_c_SEMdp);
//     omega_c_SEMdp = omega_c_SEMdp * au_2_ps * 10.0 / 53.09;
//     fscanf(idinp, "%d", &setm->N_bath_SEMdp);
    
//     setm->Nstate_SEMdp = 8;
//     *Nstate = setm->Nstate_SEMdp;
//     *Ndof1 = setm->Nstate_SEMdp - 1;
//     *Ndof2 = setm->N_bath_SEMdp;

//     H_ele = (double *)malloc(8 * 8 * sizeof(double));
//     for (int i = 0; i < 8 * 8; i++) H_ele[i] = 0.0;

//     double H_ele_values[7][7] = {
//         {12410.0, -87.7, 5.5, -5.9, 6.7, -13.7, -9.9},
//         {-87.7, 12530.0, 30.8, 8.2, 0.7, 11.8, 4.3},
//         {5.5, 30.8, 12210.0, -53.5, -2.2, -9.6, 6.0},
//         {-5.9, 8.2, -53.5, 12320.0, -70.7, -17.0, -63.3},
//         {6.7, 0.7, -2.2, -70.7, 12480.0, 81.1, -1.3},
//         {-13.7, 11.8, -9.6, -17.0, 81.1, 12630.0, 39.7},
//         {-9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 12440.0}
//     };

//     for (int i = 0; i < 7; i++) {
//         for (int j = 0; j < 7; j++) {
//             H_ele[i * 8 + j] = H_ele_values[i][j] / au_2_wn;
//         }
//     }

//     H_ele[7 * 8 + 7] += lambda_SEMdp;

//     setm->dipole_SEMdp = (double *)malloc(8 * sizeof(double));
//     setm->dipole_SEMdp[0] = 0.0; // Initialize setm->dipole_SEMdp array
// }

void parameter_FMOdp(double *mass, struct set_host *setm) {
    setm->H_ele_SEMdp = (double *)malloc(8 * 8 * sizeof(double));
    for (int i = 0; i < 8 * 8; i++) setm->H_ele_SEMdp[i] = 0.0;

    double H_ele_values[7][7] = {
        {12410.0, -87.7, 5.5, -5.9, 6.7, -13.7, -9.9},
        {-87.7, 12530.0, 30.8, 8.2, 0.7, 11.8, 4.3},
        {5.5, 30.8, 12210.0, -53.5, -2.2, -9.6, 6.0},
        {-5.9, 8.2, -53.5, 12320.0, -70.7, -17.0, -63.3},
        {6.7, 0.7, -2.2, -70.7, 12480.0, 81.1, -1.3},
        {-13.7, 11.8, -9.6, -17.0, 81.1, 12630.0, 39.7},
        {-9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 12440.0}
    };

    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 7; j++) {
            setm->H_ele_SEMdp[i * 8 + j] = H_ele_values[i][j] / au_2_wn;
        }
    }

    for (int i = 0; i < (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp; i++) {
        mass[i] = 1.0;
    }

  

    setm->c_SEMdp = (double *)malloc(setm->N_bath_SEMdp * sizeof(double));
    setm->omega_SEMdp = (double *)malloc(setm->N_bath_SEMdp * sizeof(double));
   

    for (int j = 1; j <= setm->N_bath_SEMdp; j++) {
        setm->omega_SEMdp[j-1] = setm->omega_c_SEMdp * tan(0.5 * M_PI * (1.0 - (double)j / (setm->N_bath_SEMdp + 1)));
        setm->c_SEMdp[j-1] = setm->omega_SEMdp[j-1] * sqrt(setm->lambda_SEMdp * 2.0 / (setm->N_bath_SEMdp + 1));
    }


    setm->dipole_SEMdp = (double *)malloc(setm->Nstate_SEMdp * sizeof(double));
    
}

void sample_FMOdp(double *P, double *R, double beta, struct set_host *setm) {
    sample_SEMdp(P, R, beta,setm);
}

void V_FMOdp(double *R, double *H, int forcetype, struct set_host *setm) {
    V_SEMdp(R, H, forcetype, setm);
}

void dV_FMOdp(double *R, double *dH, int forcetype, struct set_host *setm) {
    dV_SEMdp(R, dH, forcetype, setm);
}

void nucforce_FMOdp(double *R, double *nf, struct set_host *setm) {
    nucforce_SEMdp(R, nf, setm);
}

void cfweight_FMOdp(double *w0, double *wt, int idire, struct set_host *setm) {
    int i;
    double dp[8];

    switch (idire) {
        case 1:
            dp[0] = -0.741; dp[1] = -0.857; dp[2] = -0.197; dp[3] = -0.799;
            dp[4] = -0.737; dp[5] = -0.135; dp[6] = -0.495; dp[7] =  0.0;
            break;
        case 2:
            dp[0] = -0.561; dp[1] =  0.504; dp[2] =  0.957; dp[3] = -0.534;
            dp[4] =  0.656; dp[5] = -0.879; dp[6] = -0.708; dp[7] =  0.0;
            break;
        case 3:
            dp[0] = -0.3696; dp[1] = -0.107; dp[2] = -0.211; dp[3] = -0.277;
            dp[4] =  0.164; dp[5] =  0.457; dp[6] = -0.503; dp[7] =  0.0;
            break;
    }
    for (i = 0; i < setm->Nstate_SEMdp; i++){
        setm->dipole_SEMdp[i] = dp[i];
    }

    cfweight_SEMdp(w0,wt,setm);


}