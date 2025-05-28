#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif

// #define setm->omega_aso 5.0E-3
// #define setm->mass_aso 20000.0

// int type_aso;
// double setm->A_aso[9], setm->beta_aso[3], setm->D_aso[3], setm->C_aso[3], setm->R0_aso[3];
// double setm->R00_aso[9], setm->alpha_aso[9];
// double setm->Re_aso;

// void readinp_aso(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &type_aso);
//     *Ndof1 = 1;
//     *Ndof2 = 1;
//     *Nstate = 3;
// }

void parameter_aso(double *mass, struct set_host *setm);

void sample_aso(double *P, double *R, struct set_host *setm);

void V_aso(double *R, double complex *H, struct set_host *setm);

void dV_aso(double *R, double complex *dH, struct set_host *setm);

void nac_aso(double *R, double complex *nac, struct set_host *setm);