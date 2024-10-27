#ifndef MSMODEL_H
#define MSMODEL_H

#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <stdio.h>
#include "gmath.h"
#include "sbm.h"
#include "cJSON.h"
#include "msmodelio.h"
#include "def_host.h"


// extern int forcetype;
// extern char msmodelname[200];



// void readinp_msmodel(cJSON *json, int *Ndof1, int *Ndof2, int *Nstate) ;

void init_msmodel(double *mass, struct set_host *setm);

// Sample the initial conditionals for trajectories of the model
void sample_msmodel(double *P, double *R, double beta, struct set_host *setm);

// Build the diabatic potential matrix of the model
void V_msmodel(double *R, double *H, double t, struct set_host *setm);

// Build the first-order derivative matrix of the model
void dV_msmodel(double *R, double *dH, struct set_host *setm);
// Calculate the nuclear force of the model

void nucforce_msmodel(double *R, double *nf, struct set_host *setm);

void cfweight_msmodel(double *rho0, double *rhot, double beta, double *R, double *P, int icfall, struct set_host *setm);

#endif 
