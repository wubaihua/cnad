#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif

// #define setm->omega_isc 5.0E-3
// #define setm->mass_isc 20000.0

// int type_isc;
// double setm->A_isc[9], setm->beta_isc[3], setm->D_isc[3], setm->C_isc[3], setm->R0_isc[3];
// double setm->R00_isc[9], setm->alpha_isc[9];
// double setm->Re_isc;

// void readinp_isc(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &type_isc);
//     *Ndof1 = 1;
//     *Ndof2 = 1;
//     *Nstate = 3;
// }

void parameter_isc(double *mass, struct set_host *setm);

void sample_isc(double *P, double *R, struct set_host *setm);

void V_isc(double *R, double complex *H, struct set_host *setm);



void dV_isc(double *R, double complex *dH, struct set_host *setm);



void nac_isc(double *R, double complex *nac, struct set_host *setm);
