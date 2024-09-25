#include <stdio.h>
#include <math.h>

#include "msmodelio.h"
#include "def_host.h"
// #include <slave.h>
// #include <athread.h>

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

void parameter_morse3(double *mass, struct set_host *setm) ;

void sample_morse3(double *P, double *R, struct set_host *setm) ;

void V_morse3(double *R, double *H, struct set_host *setm) ;

void dV_morse3(double *R, double *dH, struct set_host *setm) ;