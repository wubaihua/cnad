
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <stdio.h>
#include "gmath.h"
#include "sbm.h"
#include "cJSON.h"



// extern int forcetype;
// extern char msmodelname[200];



// void readinp_msmodel(cJSON *json, int *Ndof1, int *Ndof2, int *Nstate) ;

void init_msmodel(double *mass);

// Sample the initial conditionals for trajectories of the model
void sample_msmodel(double *P, double *R, double beta);

// Build the diabatic potential matrix of the model
void V_msmodel(double *R, double *H, double t);

// Build the first-order derivative matrix of the model
void dV_msmodel(double *R, double *dH);
// Calculate the nuclear force of the model

void nucforce_msmodel(double *R, double *nf);




