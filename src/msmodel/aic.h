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

// int setm->N_mode_aic, setm->Nstate_aic;
// double setm->L_aic;
// double *setm->eps_aic, *setm->miu_aic, *setm->omega_aic, *setm->lambda_aic;

// void readinp_AIC(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &setm->Nstate_aic);
//     fscanf(idinp, "%d", &setm->N_mode_aic);

//     *Ndof1 = 1;
//     *Ndof2 = setm->N_mode_aic;
//     *Nstate = setm->Nstate_aic;
// }

void parameter_AIC(double *mass, struct set_host *setm);

void sample_AIC(double *P, double *R, struct set_host *setm) ;

void V_AIC(double *R, double *H, int forcetype, struct set_host *setm) ;

void dV_AIC(double *R, double *dH, int forcetype, struct set_host *setm) ;

void nucforce_AIC(double *R, double *nf, struct set_host *setm) ;

