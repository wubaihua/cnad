

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
#ifdef x86
    #include "def.h"
#endif
// #ifdef sunway
//     #include <slave.h>
//     #include <athread.h>
// #endif


// Initialize model parameters
void parameter_mole(double *mass, struct set_host *setm);

// Sample the initial conditionals for trajectories of the model
void sample_mole(double *P, double *R, double beta, struct set_host *setm);

// Build the diabatic potential matrix of the model
void V_mole(double *R, double *H, int forcetype, struct set_host *setm);

// Build the first-order derivative matrix of the model
void dV_mole(double *R, double complex *dH, int forcetype, struct set_host *setm);

// void qm_mole(double *R, struct set_host *setm, struct set_slave *sets);

// Calculate the nuclear force of the model
void nucforce_mole(double *R, double *nf, struct set_host *setm);