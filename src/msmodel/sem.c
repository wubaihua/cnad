#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"
#include <string.h>
#include <stdio.h>
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif


// int setm->N_bath_SEM, setm->Nstate_SEM;
// double setm->bias_SEM, setm->delta_SEM, setm->k_eff_SEM, setm->omega_setm->c_SEM, setm->lambda_SEM;
// double *setm->c_SEM, *setm->omega_SEM, *setm->H_ele_SEM;

// void readinp_SEM(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%lf", &setm->k_eff_SEM);
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%lf", &setm->bias_SEM);
//     setm->bias_SEM /= AU_2_EV;
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%lf", &setm->delta_SEM);
//     setm->delta_SEM /= AU_2_EV;
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%lf", &setm->omega_setm->c_SEM);
//     setm->omega_setm->c_SEM = setm->omega_setm->c_SEM * AU_2_PS * 10 / 53.09;
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%d", &setm->N_bath_SEM);

//     setm->lambda_SEM = setm->k_eff_SEM * setm->k_eff_SEM * setm->omega_setm->c_SEM / 2;

//     setm->Nstate_SEM = 2;
//     *Nstate = setm->Nstate_SEM;

//     *Ndof1 = setm->Nstate_SEM;
//     *Ndof2 = setm->N_bath_SEM;

//     setm->H_ele_SEM = (double *)malloc(4 * sizeof(double));
//     setm->H_ele_SEM[0] = setm->bias_SEM / 2;
//     setm->H_ele_SEM[1] = setm->delta_SEM;
//     setm->H_ele_SEM[2] = setm->delta_SEM;
//     setm->H_ele_SEM[3] = -setm->bias_SEM / 2;
// }

void parameter_SEM(double *mass, struct set_host *setm) {
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

void sample_SEM(double *P, double *R, double beta, struct set_host *setm) {
    double x2;
    for (int k = 0; k < setm->Nstate_SEM; k++) {
        for (int j = 0; j < setm->N_bath_SEM; j++) {
            if (beta > 99999) {
                box_muller(&P[k * setm->N_bath_SEM + j], &x2, sqrt(0.5 * hbar * setm->omega_SEM[j]), 0.0);
                box_muller(&R[k * setm->N_bath_SEM + j], &x2, sqrt(0.5 * hbar / setm->omega_SEM[j]), 0.0);
            } else if (beta > 0 && beta <= 99999) {
                box_muller(&P[k * setm->N_bath_SEM + j], &x2, sqrt(0.5 * hbar * setm->omega_SEM[j] / tanh(0.5 * beta * hbar * setm->omega_SEM[j])), 0.0);
                box_muller(&R[k * setm->N_bath_SEM + j], &x2, sqrt(0.5 * hbar / (tanh(0.5 * beta * hbar * setm->omega_SEM[j]) * setm->omega_SEM[j])), 0.0);
            } else if (beta < 0 && beta >= -99999) {
                box_muller(&P[k * setm->N_bath_SEM + j], &x2, 1.0 / sqrt(fabs(beta)), 0.0);
                box_muller(&R[k * setm->N_bath_SEM + j], &x2, 1.0 / (sqrt(fabs(beta)) * setm->omega_SEM[j]), 0.0);
            } else if (beta < -99999) {
                P[k * setm->N_bath_SEM + j] = 0;
                R[k * setm->N_bath_SEM + j] = 0;
            }
        }
    }
}

void V_SEM(double *R, double *H, int forcetype, struct set_host *setm) {
    // printf("%18.8E\n",R[0]);
    for (int i = 0; i < setm->Nstate_SEM * setm->Nstate_SEM; i++) {
        H[i] = setm->H_ele_SEM[i];
    }

    for (int i = 0; i < setm->Nstate_SEM; i++) {
        for (int j = 0; j < setm->N_bath_SEM; j++) {
            H[i * setm->Nstate_SEM + i] -= setm->c_SEM[j] * R[i * setm->N_bath_SEM + j];
        }
    }

    if (forcetype == 0) {
        double Vnuc = 0;
        for (int i = 0; i < setm->Nstate_SEM; i++) {
            for (int j = 0; j < setm->N_bath_SEM; j++) {
                Vnuc += 0.5 * setm->omega_SEM[j] * setm->omega_SEM[j] * R[i * setm->N_bath_SEM + j] * R[i * setm->N_bath_SEM + j];
            }
        }

        for (int i = 0; i < setm->Nstate_SEM; i++) {
            H[i * setm->Nstate_SEM + i] += Vnuc;
        }
    }

    // printf("%18.8E %18.8E %18.8E %18.8E\n",H[0],H[1],H[2],H[3]);
    
}

void dV_SEM(double *R, double *dH, int forcetype, struct set_host *setm) {
    for (int i = 0; i < setm->Nstate_SEM * setm->Nstate_SEM * setm->Nstate_SEM * setm->N_bath_SEM; i++) {
        dH[i] = 0;
    }

    for (int i = 0; i < setm->Nstate_SEM; i++) {
        for (int j = 0; j < setm->N_bath_SEM; j++) {
            dH[  i * setm->Nstate_SEM * setm->Nstate_SEM * setm->N_bath_SEM 
               + i * setm->Nstate_SEM * setm->N_bath_SEM + i * setm->N_bath_SEM + j] -= setm->c_SEM[j];
        }
    }

    if (forcetype == 0) {
        
        double *force = (double *)malloc(setm->Nstate_SEM * setm->N_bath_SEM * sizeof(double));
        for (int i = 0; i < setm->Nstate_SEM; i++) {
            for (int j = 0; j < setm->N_bath_SEM; j++) {
                force[i * setm->N_bath_SEM + j] = setm->omega_SEM[j] * setm->omega_SEM[j] * R[i * setm->N_bath_SEM + j];
            }
        }

        for (int i = 0; i < setm->Nstate_SEM; i++) {
            for (int j = 0; j < setm->Nstate_SEM; j++) {
                for (int k = 0; k < setm->N_bath_SEM; k++) {
                    dH[  i * setm->Nstate_SEM * setm->Nstate_SEM * setm->N_bath_SEM 
                    + i * setm->Nstate_SEM * setm->N_bath_SEM + j * setm->N_bath_SEM + k] += force[j * setm->N_bath_SEM + k];
                }
            }
        }
        free(force); 
    }
    
}

void nucforce_SEM(double *R, double *nf, struct set_host *setm) {
    for (int i = 0; i < setm->Nstate_SEM; i++) {
        for (int j = 0; j < setm->N_bath_SEM; j++) {
            nf[i * setm->N_bath_SEM + j] = setm->omega_SEM[j] * setm->omega_SEM[j] * R[i * setm->N_bath_SEM + j];
        }
    }
}

