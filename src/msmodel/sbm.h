#ifndef SBM_H
#define SBM_H

#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"  // Assuming you have a math.h header for functions like box_muller
#include <stdio.h>
// #include "cJSON.h"
#include "msmodelio.h"
#include "def_host.h"

// Spin-Boson Model parameters
// extern int N_bath_SBM, bathtype; // bathtype=1 for Ohmic; bathtype=2 for Debye
// extern double eps_SBM, delta_SBM, alpha_SBM, omega_c_SBM, lambda_SBM, s_SBM;
// extern double *c_SBM, *omega_SBM;

// Read the model type from the input file
// void readinp_SBM(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate);
// Initialize model parameters
void parameter_SBM(double *mass, struct set_host *setm);

// Sample the initial conditionals for trajectories of the model
void sample_SBM(double *P, double *R, double beta, struct set_host *setm);

// Build the diabatic potential matrix of the model
void V_SBM(double *R, double *H, int forcetype, struct set_host *setm);

// Build the first-order derivative matrix of the model
void dV_SBM(double *R, double *dH, int forcetype, struct set_host *setm) ;

// Calculate the nuclear force of the model
void nucforce_SBM(double *R, double *nf, struct set_host *setm) ;

// Compute the cfweight of the model
// void cfweight_SBM(double w0[2][2], double wt[2][2], double beta) ;



#endif 