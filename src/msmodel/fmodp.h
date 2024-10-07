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

void parameter_FMOdp(double *mass, struct set_host *setm);

void sample_FMOdp(double *P, double *R, double beta, struct set_host *setm) ;

void V_FMOdp(double *R, double *H, int forcetype, struct set_host *setm);

void dV_FMOdp(double *R, double *dH, int forcetype, struct set_host *setm);

void nucforce_FMOdp(double *R, double *nf, struct set_host *setm);

void cfweight_FMOdp(double *w0, double *wt, int idire, struct set_host *setm);