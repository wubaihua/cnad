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
#include "sem.h"

// double lambda_SEM, omega_c_SEM;
// int N_bath_SEM;

// void readinp_FMO(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%lf", &lambda_SEM);
//     lambda_SEM = lambda_SEM * AU_2_PS * 10 / 53.09;
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%lf", &omega_c_SEM);
//     omega_c_SEM = omega_c_SEM * AU_2_PS * 10 / 53.09;
//     fscanf(idinp, "%*s");
//     fscanf(idinp, "%d", &N_bath_SEM);

//     *Nstate = NSTATE_SEM;
//     *Ndof1 = NSTATE_SEM;
//     *Ndof2 = N_bath_SEM;
// }

void parameter_FMO(double *mass, struct set_host *setm);

void sample_FMO(double *P, double *R, double beta, struct set_host *setm);

void V_FMO(double *R, double *H, int forcetype, struct set_host *setm) ;

void dV_FMO(double *R, double *dH, int forcetype, struct set_host *setm) ;


void nucforce_FMO(double *R, double *nf, struct set_host *setm) ;

