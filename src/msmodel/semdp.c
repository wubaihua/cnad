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



// int setm->type_cal_SEMdp; // 1 for PT, 2 for NP
// int setm->N_bath_SEMdp, setm->Nstate_SEMdp;
// double setm->bias_SEMdp, setm->delta_SEMdp, setm->k_eff_SEMdp, setm->omega_c_SEMdp, setm->lambda_SEMdp;
// double *setm->c_SEMdp, *setm->omega_SEMdp, *setm->H_ele_SEMdp;
// double *setm->dipole_SEMdp;

// void readinp_SEMdp(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &setm->type_cal_SEMdp);
//     fscanf(idinp, "%lf", &setm->lambda_SEMdp);
//     setm->lambda_SEMdp = setm->lambda_SEMdp / au_2_wn;
//     fscanf(idinp, "%lf", &setm->bias_SEMdp);
//     setm->bias_SEMdp = setm->bias_SEMdp / au_2_wn;
//     fscanf(idinp, "%lf", &setm->delta_SEMdp);
//     setm->delta_SEMdp = setm->delta_SEMdp / au_2_wn;
//     fscanf(idinp, "%lf", &setm->omega_c_SEMdp);
//     setm->omega_c_SEMdp = setm->omega_c_SEMdp * au_2_ps * 10.0 / 53.09;
//     fscanf(idinp, "%d", &setm->N_bath_SEMdp);
    
//     setm->Nstate_SEMdp = 3;
//     *Nstate = setm->Nstate_SEMdp;
//     *Ndof1 = setm->Nstate_SEMdp - 1;
//     *Ndof2 = setm->N_bath_SEMdp;

//     setm->H_ele_SEMdp = (double *)malloc(3 * 3 * sizeof(double));
//     for (int i = 0; i < 3 * 3; i++) setm->H_ele_SEMdp[i] = 0.0;
//     setm->H_ele_SEMdp[0 * 3 + 0] = setm->bias_SEMdp + 1000.0 / au_2_wn;
//     setm->H_ele_SEMdp[0 * 3 + 1] = setm->delta_SEMdp;
//     setm->H_ele_SEMdp[1 * 3 + 1] = 1000.0 / au_2_wn;
//     setm->H_ele_SEMdp[1 * 3 + 0] = setm->delta_SEMdp;

//     if (setm->type_cal_SEMdp == 2) {
//         setm->H_ele_SEMdp[0 * 3 + 2] = 2000 * au_2_ps * 10.0 / 53.09;
//         setm->H_ele_SEMdp[2 * 3 + 0] = 2000 * au_2_ps * 10.0 / 53.09;
//         setm->H_ele_SEMdp[1 * 3 + 2] = -400 * au_2_ps * 10.0 / 53.09;
//         setm->H_ele_SEMdp[2 * 3 + 1] = -400 * au_2_ps * 10.0 / 53.09;
//     }

//     setm->dipole_SEMdp = (double *)malloc(3 * sizeof(double));
//     setm->dipole_SEMdp[0] = 5.0;
//     setm->dipole_SEMdp[1] = -1.0;
//     setm->dipole_SEMdp[2] = 0.0;
// }


void parameter_SEMdp(double *mass, struct set_host *setm) {
    int j;

    for (j = 0; j < (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp; j++) {
        mass[j] = 1.0;
    }

    setm->c_SEMdp = (double *)malloc(setm->N_bath_SEMdp * sizeof(double));
    setm->omega_SEMdp = (double *)malloc(setm->N_bath_SEMdp * sizeof(double));
    
   

    for (j = 1; j <= setm->N_bath_SEMdp; j++) {
        setm->omega_SEMdp[j-1] = setm->omega_c_SEMdp * tan(0.5 * M_PI * (1.0 - (double)j / (setm->N_bath_SEMdp + 1)));
        setm->c_SEMdp[j-1] = setm->omega_SEMdp[j-1] * sqrt(setm->lambda_SEMdp * 2.0 / (setm->N_bath_SEMdp + 1));
    }

     

    setm->H_ele_SEMdp = (double *)malloc(setm->Nstate_SEMdp * setm->Nstate_SEMdp * sizeof(double));

    for (int i = 0; i < 3 * 3; i++) setm->H_ele_SEMdp[i] = 0.0;
    setm->H_ele_SEMdp[0 * 3 + 0] = setm->bias_SEMdp + 1000.0 / au_2_wn;
    setm->H_ele_SEMdp[0 * 3 + 1] = setm->delta_SEMdp;
    setm->H_ele_SEMdp[1 * 3 + 1] = 1000.0 / au_2_wn;
    setm->H_ele_SEMdp[1 * 3 + 0] = setm->delta_SEMdp;

    if (setm->type_cal_SEMdp == 2) {
        setm->H_ele_SEMdp[0 * 3 + 2] = 2000 * au_2_ps * 10.0 / 53.09;
        setm->H_ele_SEMdp[2 * 3 + 0] = 2000 * au_2_ps * 10.0 / 53.09;
        setm->H_ele_SEMdp[1 * 3 + 2] = -400 * au_2_ps * 10.0 / 53.09;
        setm->H_ele_SEMdp[2 * 3 + 1] = -400 * au_2_ps * 10.0 / 53.09;
    }

    setm->dipole_SEMdp = (double *)malloc(3 * sizeof(double));
    setm->dipole_SEMdp[0] = 5.0;
    setm->dipole_SEMdp[1] = -1.0;
    setm->dipole_SEMdp[2] = 0.0;
    
}

void sample_SEMdp(double *P, double *R, double beta, struct set_host *setm) {
    int j, k;
    double x2;

    for (k = 0; k < setm->Nstate_SEMdp - 1; k++) {
        for (j = 0; j < setm->N_bath_SEMdp; j++) {
            if (beta > 99999) {
                box_muller(&P[k * setm->N_bath_SEMdp + j], &x2, sqrt(0.5 * hbar * setm->omega_SEMdp[j]), 0.0);
                box_muller(&R[k * setm->N_bath_SEMdp + j], &x2, sqrt(0.5 * hbar / setm->omega_SEMdp[j]), 0.0);
            } else {
                box_muller(&P[k * setm->N_bath_SEMdp + j], &x2, sqrt(0.5 * hbar * setm->omega_SEMdp[j] / tanh(0.5 * beta * hbar * setm->omega_SEMdp[j])), 0.0);
                box_muller(&R[k * setm->N_bath_SEMdp + j], &x2, sqrt(0.5 * hbar / (tanh(0.5 * beta * hbar * setm->omega_SEMdp[j]) * setm->omega_SEMdp[j])), 0.0);
                
            }
        }
    }
}

void V_SEMdp(double *R, double *H, int forcetype, struct set_host *setm) {
    int i,j;
    double Vnuc,sum;
    int n = setm->Nstate_SEMdp * setm->Nstate_SEMdp;

    for (i = 0; i < n; i++) {
        H[i] = setm->H_ele_SEMdp[i];
    }

    sum=0;
    for (i = 0; i < setm->N_bath_SEMdp; i++) {
        sum += 0.5 * setm->c_SEMdp[i] * setm->c_SEMdp[i] / (setm->omega_SEMdp[i] * setm->omega_SEMdp[i]);
    }
    for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
        H[i * setm->Nstate_SEMdp + i] += sum;
        for (j = 0; j < setm->N_bath_SEMdp; j++) {
            H[i * setm->Nstate_SEMdp + i] += -1.0 * setm->c_SEMdp[j] * R[i *  setm->N_bath_SEMdp + j];
        }
    }

    // for (i=0;i<setm->N_bath_SEMdp; i++) printf("c=%18.8E\n",setm->c_SEMdp[i]);
    // for (i=0;i<setm->N_bath_SEMdp; i++) printf("w=%18.8E\n",setm->omega_SEMdp[i]);
    // printf("x=%18.8E\n",sum);
    // printf("l=%18.8E\n",setm->lambda_SEMdp);
    // exit(-1);

    if (forcetype == 0) {
        Vnuc = 0.0;
        for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
            for (j = 0; j < setm->N_bath_SEMdp; j++) {
                Vnuc += 0.5 * setm->omega_SEMdp[j] * setm->omega_SEMdp[j] * R[i *  setm->N_bath_SEMdp + j] * R[i *  setm->N_bath_SEMdp + j];
            }
        }

        for (i = 0; i < setm->Nstate_SEMdp; i++) {
            H[i * setm->Nstate_SEMdp + i] += Vnuc;
        }
    }


    // for (i = 0; i < setm->Nstate_SEMdp * setm->Nstate_SEMdp; i++) printf("%f ",setm->H_ele_SEMdp[i]);
    
   
    
}

void dV_SEMdp(double *R, double *dH, int forcetype, struct set_host *setm) {
    int i, j, k;
    int n = setm->Nstate_SEMdp * setm->Nstate_SEMdp * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp;
    double *force = (double *)malloc((setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp * sizeof(double));

    for (i = 0; i < n; i++) {
        dH[i] = 0.0;
    }

    if (forcetype == 0) {
        for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
            for (j = 0; j < setm->N_bath_SEMdp; j++) {
                dH[i * setm->Nstate_SEMdp * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp 
                + i * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp  + i * setm->N_bath_SEMdp + j] -= setm->c_SEMdp[j];
            }
        }

        for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
            for (j = 0; j < setm->N_bath_SEMdp; j++) {
                force[i * setm->N_bath_SEMdp + j] = setm->omega_SEMdp[j] * setm->omega_SEMdp[j] * R[i * setm->N_bath_SEMdp + j];
            }
        }

        for (i = 0; i < setm->Nstate_SEMdp; i++) {
            for (j = 0; j < setm->Nstate_SEMdp - 1; j++) {
                for (k = 0; k < setm->N_bath_SEMdp; k++) {
                dH[i * setm->Nstate_SEMdp * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp 
                + i * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp + j *  setm->N_bath_SEMdp + k] += force[j *  setm->N_bath_SEMdp + k];
                }
            }
        }
        // for (i = 0; i < setm->Nstate_SEMdp; i++) {
        //     for (j = 0; j < (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp; j++) {
        //         dH[i * setm->Nstate_SEMdp * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp + j] += force[j];
        //     }
        // }
    } else if (forcetype == 1) {
        for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
            for (j = 0; j < setm->N_bath_SEMdp; j++) {
                dH[i * setm->Nstate_SEMdp * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp 
                + i * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp  + i * setm->N_bath_SEMdp + j] -= setm->c_SEMdp[j];
            }
        }
    }

    free(force);

    // for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
    //     for (j = 0; j < setm->N_bath_SEMdp; j++) {
    //         printf("%d %d %f\n",i,j,dH[i * setm->Nstate_SEMdp * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp 
    //         + i * (setm->Nstate_SEMdp - 1) * setm->N_bath_SEMdp  + i * setm->N_bath_SEMdp + j]);
    //     }
    // }
    // exit(-1);
}

void nucforce_SEMdp(double *R, double *nf, struct set_host *setm) {
    int i, j;

    for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
        for (j = 0; j < setm->N_bath_SEMdp; j++) {
            nf[i * setm->N_bath_SEMdp + j] = setm->omega_SEMdp[j] * setm->omega_SEMdp[j] * R[i * setm->N_bath_SEMdp + j];
        }
    }
}

void cfweight_SEMdp(double *w0, double *wt, struct set_host *setm) {
    int i, j;
    int n = setm->Nstate_SEMdp * setm->Nstate_SEMdp;

    for (i = 0; i < n; i++) {
        w0[i] = 0.0;
        wt[i] = 0.0;
    }

    for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
        w0[i * setm->Nstate_SEMdp + (setm->Nstate_SEMdp - 1)] = setm->dipole_SEMdp[i];
        wt[i * setm->Nstate_SEMdp + (setm->Nstate_SEMdp - 1)] = setm->dipole_SEMdp[i];
        wt[(setm->Nstate_SEMdp - 1) * setm->Nstate_SEMdp + i] = setm->dipole_SEMdp[i];
    }
}
