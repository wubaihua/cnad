
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <stdio.h>
#include "gmath.h"
#include "cJSON.h"
#include "msmodelio.h"
#include "def_host.h"

#include "sbm.h"
#include "morse3.h"
#include "sem.h"
#include "fmo.h"
#include "sf.h"
#include "aic.h"
#include "pyrazine.h"
#include "crco5.h"
#include "tully.h"
#include "semdp.h"
#include "fmodp.h"
#include "rubrene.h"
#include "dnalvcm.h"
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
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        parameter_morse3(mass,setm);
    } else if (strcmp(setm->msmodelname, "SEM") == 0 ||
       strcmp(setm->msmodelname, "sem") == 0) {
        parameter_SEM(mass,setm);
    } else if (strcmp(setm->msmodelname, "FMO") == 0 ||
       strcmp(setm->msmodelname, "fmo") == 0) {
        parameter_FMO(mass,setm);
    } else if (strcmp(setm->msmodelname, "SF") == 0 ||
       strcmp(setm->msmodelname, "sf") == 0) {
        parameter_SF(mass,setm);
    } else if (strcmp(setm->msmodelname, "AIC") == 0 ||
       strcmp(setm->msmodelname, "aic") == 0) {
        parameter_AIC(mass,setm);
    } else if (strcmp(setm->msmodelname, "pyrazine") == 0) {
        parameter_pyrazine(mass,setm);
    } else if (strcmp(setm->msmodelname, "crco5") == 0 ||
       strcmp(setm->msmodelname, "CrCO5") == 0) {
        parameter_crco5(mass,setm);
    } else if (strcmp(setm->msmodelname, "tully") == 0) {
        parameter_tully(mass,setm);
    } else if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        parameter_SEMdp(mass,setm);
    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        parameter_FMOdp(mass,setm);
    } else if (strcmp(setm->msmodelname, "rubrene") == 0) {
        parameter_rubrene(mass,setm);
    } else if (strcmp(setm->msmodelname, "dnalvcm") == 0) {
        parameter_dnalvcm(mass,setm);
    }

}

// Sample the initial conditionals for trajectories of the model
void sample_msmodel(double *P, double *R, double beta, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        sample_SBM(P, R, beta,setm);
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        sample_morse3(P, R,setm);
    } else if (strcmp(setm->msmodelname, "SEM") == 0 ||
       strcmp(setm->msmodelname, "sem") == 0) {
        sample_SEM(P, R, beta,setm);
    } else if (strcmp(setm->msmodelname, "FMO") == 0 ||
       strcmp(setm->msmodelname, "fmo") == 0) {
        sample_FMO(P, R, beta,setm);
    } else if (strcmp(setm->msmodelname, "SF") == 0 ||
       strcmp(setm->msmodelname, "sf") == 0) {
        sample_SF(P, R, beta,setm);
    } else if (strcmp(setm->msmodelname, "AIC") == 0 ||
       strcmp(setm->msmodelname, "aic") == 0) {
        sample_AIC(P, R, setm);
    } else if (strcmp(setm->msmodelname, "pyrazine") == 0) {
        sample_pyrazine(P, R, setm);
    } else if (strcmp(setm->msmodelname, "crco5") == 0 ||
       strcmp(setm->msmodelname, "CrCO5") == 0) {
        sample_crco5(P, R, setm);
    } else if (strcmp(setm->msmodelname, "tully") == 0) {
        sample_tully(P, R, setm);
    } else if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        sample_SEMdp(P, R, beta,setm);
    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        sample_FMOdp(P, R, beta,setm);
    } else if (strcmp(setm->msmodelname, "rubrene") == 0) {
        sample_rubrene(P, R, beta,setm);
    } else if (strcmp(setm->msmodelname, "dnalvcm") == 0) {
        sample_dnalvcm(P, R, setm);
    }
}

// Build the diabatic potential matrix of the model
void V_msmodel(double *R, double *H, double t, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        V_SBM(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        V_morse3(R, H, setm);
    } else if (strcmp(setm->msmodelname, "SEM") == 0 ||
       strcmp(setm->msmodelname, "sem") == 0) {
        V_SEM(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "FMO") == 0 ||
       strcmp(setm->msmodelname, "fmo") == 0) {
        V_FMO(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "SF") == 0 ||
       strcmp(setm->msmodelname, "sf") == 0) {
        V_SF(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "AIC") == 0 ||
       strcmp(setm->msmodelname, "aic") == 0) {
        V_AIC(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "pyrazine") == 0) {
        V_pyrazine(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "crco5") == 0 ||
       strcmp(setm->msmodelname, "CrCO5") == 0) {
        V_crco5(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "tully") == 0) {
        V_tully(R, H, setm);
    } else if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        V_SEMdp(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        V_FMOdp(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "rubrene") == 0) {
        V_rubrene(R, H, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "dnalvcm") == 0 ) {
        V_dnalvcm(R, H, setm->forcetype,setm);
    }
}

// Build the first-order derivative matrix of the model
void dV_msmodel(double *R, double *dH, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        dV_SBM(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        dV_morse3(R, dH, setm);
    } else if (strcmp(setm->msmodelname, "SEM") == 0 ||
       strcmp(setm->msmodelname, "sem") == 0) {
        dV_SEM(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "FMO") == 0 ||
       strcmp(setm->msmodelname, "fmo") == 0) {
        dV_FMO(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "SF") == 0 ||
       strcmp(setm->msmodelname, "sf") == 0) {
        dV_SF(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "AIC") == 0 ||
       strcmp(setm->msmodelname, "aic") == 0) {
        dV_AIC(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "pyrazine") == 0) {
        dV_pyrazine(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "crco5") == 0 ||
       strcmp(setm->msmodelname, "CrCO5") == 0) {
        dV_crco5(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "tully") == 0) {
        dV_tully(R, dH, setm);
    } else if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        dV_SEMdp(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        dV_FMOdp(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "rubrene") == 0) {
        dV_rubrene(R, dH, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "dnalvcm") == 0 ) {
        dV_dnalvcm(R, dH, setm->forcetype,setm);
    }
}

// Calculate the nuclear force of the model
void nucforce_msmodel(double *R, double *nf, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        nucforce_SBM(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "SEM") == 0 ||
       strcmp(setm->msmodelname, "sem") == 0) {
        nucforce_SEM(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "FMO") == 0 ||
       strcmp(setm->msmodelname, "fmo") == 0) {
        nucforce_FMO(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "SF") == 0 ||
       strcmp(setm->msmodelname, "sf") == 0) {
        nucforce_SF(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "AIC") == 0 ||
       strcmp(setm->msmodelname, "aic") == 0) {
        nucforce_AIC(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "pyrazine") == 0) {
        nucforce_pyrazine(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "crco5") == 0 ||
       strcmp(setm->msmodelname, "CrCO5") == 0) {
        nucforce_crco5(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        nucforce_SEMdp(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        nucforce_FMOdp(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "rubrene") == 0) {
        nucforce_rubrene(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "dnalvcm") == 0) {
        nucforce_dnalvcm(R, nf,setm);
    }
}


void cfweight_msmodel(double *rho0, double *rhot, double beta, double *R, double *P, int icfall, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        cfweight_SEMdp(rho0,rhot,setm);
    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        cfweight_FMOdp(rho0,rhot, icfall, setm);
    } else if (strcmp(setm->msmodelname, "rubrene") == 0) {
        cfweight_rubrene(rho0,rhot, beta, R, P, setm);
    }
}






