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

// Constants and parameters
// int setm->Nstate_rubrene, setm->N_mode_rubrene, setm->Nmole_rubrene;
// double setm->Vc_rubrene, setm->factor_freq_rubrene;
// double *setm->omega_rubrene, *setm->lambda_rubrene, *setm->g_rubrene, *setm->mass_rubrene;

// void readinp_rubrene(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     char c50[50];
    
//     fscanf(idinp, "%s", c50); // Read character string
//     fscanf(idinp, "%d", &setm->Nmole_rubrene);
//     fscanf(idinp, "%s", c50);
//     fscanf(idinp, "%lf", &setm->Vc_rubrene);
//     fscanf(idinp, "%s", c50);
//     fscanf(idinp, "%lf", &setm->factor_freq_rubrene);

//     setm->Vc_rubrene = setm->Vc_rubrene / au_2_eV;

//     setm->Nstate_rubrene = setm->Nmole_rubrene;
//     setm->N_mode_rubrene = 9;

//     *Ndof1 = setm->Nmole_rubrene;
//     *Ndof2 = setm->N_mode_rubrene;
//     *Nstate = setm->Nstate_rubrene;
// }




void parameter_rubrene(double *mass, struct set_host *setm);
void sample_rubrene(double *P, double *R, double beta, struct set_host *setm);

void V_rubrene(double *R, double *H, int forcetype, struct set_host *setm);

void dV_rubrene(double *R, double *dH, int forcetype, struct set_host *setm);

void nucforce_rubrene(double *R, double *nf, struct set_host *setm);




void cfweight_rubrene(double *w0, double *wt, double beta, double *R, double *P, struct set_host *setm);


