
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



// int forcetype;
// char setm->msmodelname[200];



// void readinp_msmodel(cJSON *json, int *Ndof1, int *Ndof2, int *Nstate) {
//     if (strcmp(setm->msmodelname, "SBM") == 0 ||
//        strcmp(setm->msmodelname, "sbm") == 0) {
//         readinp_SBM(json, Ndof1, Ndof2, Nstate);
//     }
// }




void init_msmodel(double *mass, struct set_host *setm){
    
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        parameter_SBM(mass,setm);
    }
}

// Sample the initial conditionals for trajectories of the model
void sample_msmodel(double *P, double *R, double beta, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        sample_SBM(P, R, beta,setm);
    }
}

// Build the diabatic potential matrix of the model
void V_msmodel(double *R, double *H, double t, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        V_SBM(R, H, setm->forcetype,setm);
    }
}

// Build the first-order derivative matrix of the model
void dV_msmodel(double *R, double *dH, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        dV_SBM(R, dH, setm->forcetype,setm);
    }
}

// Calculate the nuclear force of the model
void nucforce_msmodel(double *R, double *nf, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        nucforce_SBM(R, nf,setm);
        
    }
}






