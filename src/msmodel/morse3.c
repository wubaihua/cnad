#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#include <slave.h>
#include <athread.h>

// #define setm->omega_morse3 5.0E-3
// #define setm->mass_morse3 20000.0

// int type_morse3;
// double setm->A_morse3[9], setm->beta_morse3[3], setm->D_morse3[3], setm->C_morse3[3], setm->R0_morse3[3];
// double setm->R00_morse3[9], setm->alpha_morse3[9];
// double setm->Re_morse3;

// void readinp_morse3(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &type_morse3);
//     *Ndof1 = 1;
//     *Ndof2 = 1;
//     *Nstate = 3;
// }

void parameter_morse3(double *mass, struct set_host *setm) {
    int i, j;

    setm->mass_morse3 = 20000;
    setm->omega_morse3 = 5.0E-3;

    mass[0] = setm->mass_morse3;

    switch (setm->type_morse3) {
        case 1:
            setm->D_morse3[0] = 0.003; setm->D_morse3[1] = 0.004; setm->D_morse3[2] = 0.003;
            setm->beta_morse3[0] = 0.65; setm->beta_morse3[1] = 0.6; setm->beta_morse3[2] = 0.65;
            setm->R0_morse3[0] = 5.0; setm->R0_morse3[1] = 4.0; setm->R0_morse3[2] = 6.0;
            setm->C_morse3[0] = 0.0; setm->C_morse3[1] = 0.01; setm->C_morse3[2] = 0.006;
            setm->A_morse3[1] = 0.002; setm->A_morse3[5] = 0.002; setm->A_morse3[2] = 0;
            setm->R00_morse3[1] = 3.4; setm->R00_morse3[5] = 4.8; setm->R00_morse3[2] = 0;
            setm->alpha_morse3[1] = 16; setm->alpha_morse3[5] = 16; setm->alpha_morse3[2] = 0;
            setm->Re_morse3 = 2.9;
            break;
        case 2:
            setm->D_morse3[0] = 0.02; setm->D_morse3[1] = 0.01; setm->D_morse3[2] = 0.003;
            setm->beta_morse3[0] = 0.65; setm->beta_morse3[1] = 0.4; setm->beta_morse3[2] = 0.65;
            setm->R0_morse3[0] = 4.5; setm->R0_morse3[1] = 4.0; setm->R0_morse3[2] = 4.4;
            setm->C_morse3[0] = 0.0; setm->C_morse3[1] = 0.01; setm->C_morse3[2] = 0.02;
            setm->A_morse3[1] = 0.005; setm->A_morse3[5] = 0; setm->A_morse3[2] = 0.005;
            setm->R00_morse3[1] = 3.66; setm->R00_morse3[5] = 0; setm->R00_morse3[2] = 3.34;
            setm->alpha_morse3[1] = 32; setm->alpha_morse3[5] = 0; setm->alpha_morse3[2] = 32;
            setm->Re_morse3 = 3.3;
            break;
        case 3:
            setm->D_morse3[0] = 0.02; setm->D_morse3[1] = 0.02; setm->D_morse3[2] = 0.003;
            setm->beta_morse3[0] = 0.4; setm->beta_morse3[1] = 0.65; setm->beta_morse3[2] = 0.65;
            setm->R0_morse3[0] = 4.0; setm->R0_morse3[1] = 4.5; setm->R0_morse3[2] = 6.0;
            setm->C_morse3[0] = 0.02; setm->C_morse3[1] = 0.0; setm->C_morse3[2] = 0.02;
            setm->A_morse3[1] = 0.005; setm->A_morse3[5] = 0; setm->A_morse3[2] = 0.005;
            setm->R00_morse3[1] = 3.4; setm->R00_morse3[5] = 0; setm->R00_morse3[2] = 4.97;
            setm->alpha_morse3[1] = 32; setm->alpha_morse3[5] = 0; setm->alpha_morse3[2] = 32;
            setm->Re_morse3 = 2.1;
            break;
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            if (i > j) {
                setm->A_morse3[i * 3 + j] = setm->A_morse3[j * 3 + i];
                setm->R00_morse3[i * 3 + j] = setm->R00_morse3[j * 3 + i];
                setm->alpha_morse3[i * 3 + j] = setm->alpha_morse3[j * 3 + i];
            }
        }
    }

}

void sample_morse3(double *P, double *R, struct set_host *setm) {
    double x2, gamma_morse3;
    gamma_morse3 = setm->mass_morse3 * setm->omega_morse3;

    // Assuming box_muller is defined elsewhere
    box_muller(&P[0], &x2, sqrt(gamma_morse3 / 2), 0.0);
    box_muller(&R[0], &x2, sqrt(1 / (gamma_morse3 * 2)), setm->Re_morse3);
}

void V_morse3(double *R, double *H, struct set_host *setm) {
    int i, j;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            if (i == j) {
                H[i * 3 + i] = setm->D_morse3[i] * pow(1 - exp(-setm->beta_morse3[i] * (R[0] - setm->R0_morse3[i])), 2) + setm->C_morse3[i];
            } else {
                H[i * 3 + j] = setm->A_morse3[i * 3 + j] * exp(-setm->alpha_morse3[i * 3 + j] * pow(R[0] - setm->R00_morse3[i * 3 + j], 2));
            }
        }
    }

    // printf("%18.8E %18.8E\n",setm->beta_morse3[0], setm->R0_morse3[0]);
    // printf("%18.8E %18.8E %18.8E %18.8E\n",R[0],H[0], H[4], H[8]);

}

void dV_morse3(double *R, double *dH, struct set_host *setm) {
    int i, j;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            if (i == j) {
                dH[i * 3 + i] = 2 * setm->beta_morse3[i] * setm->D_morse3[i] * exp(-setm->beta_morse3[i] * (R[0] - setm->R0_morse3[i])) * (1 - exp(-setm->beta_morse3[i] * (R[0] - setm->R0_morse3[i])));
            } else {
                dH[i * 3 + j] = -2 * setm->A_morse3[i * 3 + j] * setm->alpha_morse3[i * 3 + j] * (R[0] - setm->R00_morse3[i * 3 + j]) * exp(-setm->alpha_morse3[i * 3 + j] * pow(R[0] - setm->R00_morse3[i * 3 + j], 2));
            }
        }
    }

}
