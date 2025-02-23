#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif

// #define setm->omega_morse2 5.0E-3
// #define setm->mass_morse2 20000.0

// int type_morse2;
// double setm->A_morse2[9], setm->beta_morse2[3], setm->D_morse2[3], setm->C_morse2[3], setm->R0_morse2[3];
// double setm->R00_morse2[9], setm->alpha_morse2[9];
// double setm->Re_morse2;

// void readinp_morse2(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &type_morse2);
//     *Ndof1 = 1;
//     *Ndof2 = 1;
//     *Nstate = 3;
// }

void parameter_morse2(double *mass, struct set_host *setm) {
    int i, j;

    setm->mass_morse2 = 20000;
    setm->omega_morse2 = 5.0E-3;

    mass[0] = setm->mass_morse2;

    switch (setm->type_morse2) {
        case 1:
            // setm->D_morse2[0] = 0.003; setm->D_morse2[1] = 0.004; setm->D_morse2[2] = 0.003;
            // setm->beta_morse2[0] = 0.65; setm->beta_morse2[1] = 0.6; setm->beta_morse2[2] = 0.65;
            // setm->R0_morse2[0] = 5.0; setm->R0_morse2[1] = 4.0; setm->R0_morse2[2] = 6.0;
            // setm->C_morse2[0] = 0.0; setm->C_morse2[1] = 0.01; setm->C_morse2[2] = 0.006;
            // setm->A_morse2[1] = 0.002; setm->A_morse2[5] = 0.002; setm->A_morse2[2] = 0;
            // setm->R00_morse2[1] = 3.4; setm->R00_morse2[5] = 4.8; setm->R00_morse2[2] = 0;
            // setm->alpha_morse2[1] = 16; setm->alpha_morse2[5] = 16; setm->alpha_morse2[2] = 0;
            // setm->Re_morse2 = 2.9;
            break;
        case 2:
            setm->D_morse2[0] = 0.02; setm->D_morse2[1] = 0.003;
            setm->beta_morse2[0] = 0.65; setm->beta_morse2[1] = 0.65;
            setm->R0_morse2[0] = 4.5; setm->R0_morse2[1] = 4.4;
            setm->C_morse2[0] = 0.0; setm->C_morse2[1] = 0.02;
            setm->A_morse2[1] = 0.005; setm->A_morse2[2] = 0.005;
            setm->R00_morse2[1] = 3.34; setm->R00_morse2[2] = 3.34;
            setm->alpha_morse2[1] = 32; setm->alpha_morse2[2] = 32;
            setm->Re_morse2 = 3.3;
            break;
        case 3:
            // setm->D_morse2[0] = 0.02; setm->D_morse2[1] = 0.02; setm->D_morse2[2] = 0.003;
            // setm->beta_morse2[0] = 0.4; setm->beta_morse2[1] = 0.65; setm->beta_morse2[2] = 0.65;
            // setm->R0_morse2[0] = 4.0; setm->R0_morse2[1] = 4.5; setm->R0_morse2[2] = 6.0;
            // setm->C_morse2[0] = 0.02; setm->C_morse2[1] = 0.0; setm->C_morse2[2] = 0.02;
            // setm->A_morse2[1] = 0.005; setm->A_morse2[5] = 0; setm->A_morse2[2] = 0.005;
            // setm->R00_morse2[1] = 3.4; setm->R00_morse2[5] = 0; setm->R00_morse2[2] = 4.97;
            // setm->alpha_morse2[1] = 32; setm->alpha_morse2[5] = 0; setm->alpha_morse2[2] = 32;
            // setm->Re_morse2 = 2.1;
            break;
    }

    // for (i = 0; i < 3; i++) {
    //     for (j = 0; j < 3; j++) {
    //         if (i > j) {
    //             setm->A_morse2[i * 3 + j] = setm->A_morse2[j * 3 + i];
    //             setm->R00_morse2[i * 3 + j] = setm->R00_morse2[j * 3 + i];
    //             setm->alpha_morse2[i * 3 + j] = setm->alpha_morse2[j * 3 + i];
    //         }
    //     }
    // }

}

void sample_morse2(double *P, double *R, struct set_host *setm) {
    double x2, gamma_morse2;
    gamma_morse2 = setm->mass_morse2 * setm->omega_morse2;

    // Assuming box_muller is defined elsewhere
    box_muller(&P[0], &x2, sqrt(gamma_morse2 / 2), 0.0);
    box_muller(&R[0], &x2, sqrt(1 / (gamma_morse2 * 2)), setm->Re_morse2);
}

void V_morse2(double *R, double *H, struct set_host *setm) {
    int i, j;

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            if (i == j) {
                H[i * 2 + i] = setm->D_morse2[i] * pow(1 - exp(-setm->beta_morse2[i] * (R[0] - setm->R0_morse2[i])), 2) + setm->C_morse2[i];
            } else {
                H[i * 2 + j] = setm->A_morse2[i * 2 + j] * exp(-setm->alpha_morse2[i * 2 + j] * pow(R[0] - setm->R00_morse2[i * 2 + j], 2));
            }
        }
    }

    // printf("%18.8E %18.8E\n",setm->beta_morse2[0], setm->R0_morse2[0]);
    // printf("%18.8E %18.8E %18.8E %18.8E\n",R[0],H[0], H[4], H[8]);

}

void dV_morse2(double *R, double *dH, struct set_host *setm) {
    int i, j;

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            if (i == j) {
                dH[i * 2 + i] = 2 * setm->beta_morse2[i] * setm->D_morse2[i] * exp(-setm->beta_morse2[i] * (R[0] - setm->R0_morse2[i])) * (1 - exp(-setm->beta_morse2[i] * (R[0] - setm->R0_morse2[i])));
            } else {
                dH[i * 2 + j] = -2 * setm->A_morse2[i * 2 + j] * setm->alpha_morse2[i * 2 + j] * (R[0] - setm->R00_morse2[i * 2 + j]) * exp(-setm->alpha_morse2[i * 2 + j] * pow(R[0] - setm->R00_morse2[i * 2 + j], 2));
            }
        }
    }

}
