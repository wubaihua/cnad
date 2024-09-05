
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <stdio.h>
#include "gmath.h"
#include "sbm.h"



int forcetype;
char msmodelname[200];



void readinp_msmodel(FILE *idinp, int Ndof1, int Ndof2, int Nstate) {
    if (strcmp(trim(adjustl(msmodelname)), "SBM") == 0 ||
       strcmp(trim(adjustl(msmodelname)), "sbm") == 0) {
        readinp_SBM(idinp, Ndof1, Ndof2, Nstate);
    }
}




void parameter_msmodel(double *mass){
    if (strcmp(trim(adjustl(msmodelname)), "SBM") == 0 ||
       strcmp(trim(adjustl(msmodelname)), "sbm") == 0) {
        parameter_msmodel(mass);
    }
}

// Sample the initial conditionals for trajectories of the model
void sample_msmodel(double *P, double *R, double beta){
    if (strcmp(trim(adjustl(msmodelname)), "SBM") == 0 ||
       strcmp(trim(adjustl(msmodelname)), "sbm") == 0) {
        sample_SBM(P, R, beta);
    }
}

// Build the diabatic potential matrix of the model
void V_msmodel(double *R, double *H, int forcetype){
    if (strcmp(trim(adjustl(msmodelname)), "SBM") == 0 ||
       strcmp(trim(adjustl(msmodelname)), "sbm") == 0) {
        V_SBM(R, H, forcetype);
    }
}

// Build the first-order derivative matrix of the model
void dV_msmodel(double *R, double *dH, int forcetype){
    if (strcmp(trim(adjustl(msmodelname)), "SBM") == 0 ||
       strcmp(trim(adjustl(msmodelname)), "sbm") == 0) {
        dV_SBM(R, dH, forcetype);
    }
}

// Calculate the nuclear force of the model
void nucforce_msmodel(double *R, double *nf){
    if (strcmp(trim(adjustl(msmodelname)), "SBM") == 0 ||
       strcmp(trim(adjustl(msmodelname)), "sbm") == 0) {
        nucforce_SBM(R, nf);
    }
}






