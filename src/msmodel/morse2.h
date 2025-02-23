#include <stdio.h>
#include <math.h>

#include "msmodelio.h"
#include "def_host.h"
// #include <slave.h>
// #include <athread.h>

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

void parameter_morse2(double *mass, struct set_host *setm) ;

void sample_morse2(double *P, double *R, struct set_host *setm) ;

void V_morse2(double *R, double *H, struct set_host *setm) ;

void dV_morse2(double *R, double *dH, struct set_host *setm) ;