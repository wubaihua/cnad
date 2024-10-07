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
// double setm->bias_SEMdp, setm->delta_SEMdp, setm->k_eff_SEMdp, setm->omega_setm->c_SEMdp, setm->lambda_SEMdp;
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
//     fscanf(idinp, "%lf", &setm->omega_setm->c_SEMdp);
//     setm->omega_setm->c_SEMdp = setm->omega_setm->c_SEMdp * au_2_ps * 10.0 / 53.09;
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


void parameter_SEMdp(double *mass, struct set_host *setm) ;

void sample_SEMdp(double *P, double *R, double beta, struct set_host *setm) ;

void V_SEMdp(double *R, double *H, int forcetype, struct set_host *setm);

void dV_SEMdp(double *R, double *dH, int forcetype, struct set_host *setm) ;

void nucforce_SEMdp(double *R, double *nf, struct set_host *setm) ;

void cfweight_SEMdp(double *w0, double *wt, struct set_host *setm);