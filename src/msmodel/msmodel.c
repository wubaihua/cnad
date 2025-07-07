
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
#include "morse2.h"
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
#include "frozen.h"
#include "dualho.h"
#include "bpy.h"
#include "retinal.h"
#include "aso.h"

#ifdef x86
    #include "def.h"
    #include "mole.h"
    void qm_mole(double *R, struct set_host *setm, struct set_slave *sets);
#endif
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
    } else if (strcmp(setm->msmodelname, "morse2") == 0 ||
       strcmp(setm->msmodelname, "Morse2") == 0) {
        parameter_morse2(mass,setm);
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
    } else if (strcmp(setm->msmodelname, "frozen") == 0) {
        parameter_frozen(mass,setm);
    } else if (strcmp(setm->msmodelname, "dualho") == 0) {
        parameter_dualho(mass,setm);
    } else if (strcmp(setm->msmodelname, "bpy") == 0) {
        parameter_bpy(mass,setm);
    } else if (strcmp(setm->msmodelname, "retinal") == 0) {
        parameter_retinal(mass,setm);
    } else if (strcmp(setm->msmodelname, "aso") == 0) {
        parameter_aso(mass,setm);
    }

    #ifdef x86
        if (strcmp(setm->msmodelname, "mole") == 0) {
            parameter_mole(mass,setm);
        }
    #endif
}

// Sample the initial conditionals for trajectories of the model
void sample_msmodel(double *P, double *R, double beta, struct set_host *setm){
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        sample_SBM(P, R, beta,setm);
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        sample_morse3(P, R,setm);
    } else if (strcmp(setm->msmodelname, "morse2") == 0 ||
       strcmp(setm->msmodelname, "Morse2") == 0) {
        sample_morse2(P, R,setm);
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
    } else if (strcmp(setm->msmodelname, "frozen") == 0) {
        sample_frozen(P, R, setm);
    } else if (strcmp(setm->msmodelname, "dualho") == 0) {
        sample_dualho(P, R, setm);
    } else if (strcmp(setm->msmodelname, "bpy") == 0) {
        sample_bpy(P, R, setm);
    } else if (strcmp(setm->msmodelname, "retinal") == 0) {
        sample_retinal(P, R, setm);
    } else if (strcmp(setm->msmodelname, "aso") == 0) {
        sample_aso(P, R, setm);
    } 
    
    #ifdef x86
        if (strcmp(setm->msmodelname, "mole") == 0) {
            sample_mole(P, R, beta, setm);
        }
    #endif
    
}

// Build the diabatic potential matrix of the model
void V_msmodel(double *R, double complex *H, double t, struct set_host *setm){
    double *V_real = (double *)malloc(setm->Nstate * setm->Nstate * sizeof(double));
    int ifcpy = 0;
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        V_SBM(R, V_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        V_morse3(R, V_real, setm);
    } else if (strcmp(setm->msmodelname, "morse2") == 0 ||
       strcmp(setm->msmodelname, "Morse2") == 0) {
        V_morse2(R, V_real, setm);
    } else if (strcmp(setm->msmodelname, "SEM") == 0 ||
       strcmp(setm->msmodelname, "sem") == 0) {
        V_SEM(R, V_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "FMO") == 0 ||
       strcmp(setm->msmodelname, "fmo") == 0) {
        V_FMO(R, V_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "SF") == 0 ||
       strcmp(setm->msmodelname, "sf") == 0) {
        V_SF(R, V_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "AIC") == 0 ||
       strcmp(setm->msmodelname, "aic") == 0) {
        V_AIC(R, V_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "pyrazine") == 0) {
        V_pyrazine(R, H, setm->forcetype,setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "crco5") == 0 ||
       strcmp(setm->msmodelname, "CrCO5") == 0) {
        V_crco5(R, H, setm->forcetype,setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "tully") == 0) {
        V_tully(R, V_real, setm);
    } else if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        V_SEMdp(R, V_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        V_FMOdp(R, V_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "rubrene") == 0) {
        V_rubrene(R, V_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "dnalvcm") == 0 ) {
        V_dnalvcm(R, H, setm->forcetype,setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "frozen") == 0 ) {
        V_frozen(R, V_real, setm);
    } else if (strcmp(setm->msmodelname, "dualho") == 0 ) {
        V_dualho(R, V_real, setm);
    } else if (strcmp(setm->msmodelname, "bpy") == 0 ) {
        V_bpy(R, H, setm->forcetype,setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "retinal") == 0 ) {
        V_retinal(R, H, setm->forcetype, setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "aso") == 0 ) {
        V_aso(R, H, setm);
        ifcpy = 1;
    }

    if(ifcpy == 0){
        for (int i = 0; i < setm->Nstate; i++) {
            for (int j = 0; j < setm->Nstate; j++) {
                H[i * setm->Nstate + j] = V_real[i * setm->Nstate + j] + 0 * I;
            }
        }
    }
    free(V_real);
}

// Build the first-order derivative matrix of the model
void dV_msmodel(double *R, double complex *dH, struct set_host *setm){
    double *dV_real = (double *)malloc(setm->Nstate * setm->Nstate * setm->Ndof1 * setm->Ndof2 * sizeof(double));
    int ifcpy = 0;
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        dV_SBM(R, dV_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        dV_morse3(R, dV_real, setm);
    } else if (strcmp(setm->msmodelname, "morse2") == 0 ||
       strcmp(setm->msmodelname, "Morse2") == 0) {
        dV_morse2(R, dV_real, setm);
    } else if (strcmp(setm->msmodelname, "SEM") == 0 ||
       strcmp(setm->msmodelname, "sem") == 0) {
        dV_SEM(R, dV_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "FMO") == 0 ||
       strcmp(setm->msmodelname, "fmo") == 0) {
        dV_FMO(R, dV_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "SF") == 0 ||
       strcmp(setm->msmodelname, "sf") == 0) {
        dV_SF(R, dV_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "AIC") == 0 ||
       strcmp(setm->msmodelname, "aic") == 0) {
        dV_AIC(R, dV_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "pyrazine") == 0) {
        dV_pyrazine(R, dH, setm->forcetype,setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "crco5") == 0 ||
       strcmp(setm->msmodelname, "CrCO5") == 0) {
        dV_crco5(R, dH, setm->forcetype,setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "tully") == 0) {
        dV_tully(R, dV_real, setm);
    } else if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        dV_SEMdp(R, dV_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        dV_FMOdp(R, dV_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "rubrene") == 0) {
        dV_rubrene(R, dV_real, setm->forcetype,setm);
    } else if (strcmp(setm->msmodelname, "dnalvcm") == 0 ) {
        dV_dnalvcm(R, dH, setm->forcetype,setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "frozen") == 0 ) {
        dV_frozen(R, dV_real, setm);
    } else if (strcmp(setm->msmodelname, "dualho") == 0 ) {
        dV_dualho(R, dV_real, setm);
    } else if (strcmp(setm->msmodelname, "bpy") == 0 ) {
        dV_bpy(R, dH, setm->forcetype,setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "retinal") == 0 ) {
        dV_retinal(R, dH,setm->forcetype, setm);
        ifcpy = 1;
    } else if (strcmp(setm->msmodelname, "aso") == 0 ) {
        dV_aso(R, dH, setm);
        ifcpy = 1;
    } 
    
    // #ifdef x86
    //     if (strcmp(setm->msmodelname, "mole") == 0 ) {
    //         dV_mole(R, dH, setm->forcetype, setm);
    //         ifcpy = 1;
    //     }
    // #endif 


    if(ifcpy == 0){
        for (int i = 0; i < setm->Nstate; i++) {
            for (int j = 0; j < setm->Nstate; j++) {
                for (int k = 0; k < setm->Ndof1 * setm->Ndof2; k++) {
                    dH[i * setm->Nstate * setm->Ndof1 * setm->Ndof2 + j * setm->Ndof1 * setm->Ndof2 + k] = dV_real[i * setm->Nstate * setm->Ndof1 * setm->Ndof2 + j * setm->Ndof1 * setm->Ndof2 + k] + 0 * I;
                }
            }
        }
    }
    free(dV_real);
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
    } else if (strcmp(setm->msmodelname, "bpy") == 0) {
        nucforce_bpy(R, nf,setm);
    } else if (strcmp(setm->msmodelname, "retinal") == 0) {
        nucforce_retinal(R, nf,setm);
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
    } else if (strcmp(setm->msmodelname, "frozen") == 0) {
        cfweight_frozen(rho0,rhot,setm);
    }
}



void nac_msmodel(double *R, double complex *nac, struct set_host *setm){
    // double *dV_real = (double *)malloc(setm->Nstate * setm->Nstate * setm->Ndof1 * setm->Ndof2 * sizeof(double));
    // int ifcpy = 0;
    if (strcmp(setm->msmodelname, "aso") == 0) {
        nac_aso(R, nac, setm);
    }
}


#ifdef x86
void qm_msmodel(double *R, struct set_host *setm, struct set_slave *sets){
    // double *dV_real = (double *)malloc(setm->Nstate * setm->Nstate * setm->Ndof1 * setm->Ndof2 * sizeof(double));
    // int ifcpy = 0;
    double complex tempcm1[setm->Nstate * setm->Nstate], tempcm2[setm->Nstate * setm->Nstate];
    double tempdm1[setm->Nstate*setm->Nstate],tempdm2[setm->Nstate*setm->Nstate],tempdm3[setm->Nstate*setm->Nstate],tempdm4[setm->Nstate*setm->Nstate], tempdv1[setm->Nstate];
    int id_max[setm->Nstate], idloc, id1, id2;
    double complex overlap[setm->Nstate * setm->Nstate];
    double overlap2[setm->Nstate * setm->Nstate], vmax;
    
        if (strcmp(setm->msmodelname, "mole") == 0 ) {
            qm_mole(R, setm, sets);
        }

        // printf("V before corr:\n");
        // for (int i = 0; i < setm->Nstate; i++) {
        //     for (int j = 0; j < setm->Nstate; j++) {
        //         printf("%18.8E  ", creal(sets->V[i * setm->Nstate + j]));
        //     }
        //     printf("\n");
        // }



        
        



        if (setm->rep == 2 || setm->rep == 3) {
            for (int i = 0; i < setm->Nstate; i++) {
                for (int j = 0; j < setm->Nstate; j++) {
                    if (i == j) {
                        continue;
                    }
                    if (cabs(sets->V[i * setm->Nstate + j]) > 1e-10) {
                        if (creal(sets->V[i * setm->Nstate + j]) * (creal(sets->V_old[i * setm->Nstate + j])) < 0.0) {
                            sets->V[i * setm->Nstate + j] = - sets->V[i * setm->Nstate + j];
                        }
                    }
                }
            }
            // printf("11111\n");
            dia_hermitemat(setm->Nstate, sets->V, sets->E_adia, sets->U_d2a);
            
            if (sets->if_ad_nac) {
                // transpose(sets->U_d2a,tempdm1,setm->Nstate);
                // dd_matmul(tempdm1,sets->U_ref,overlap,setm->Nstate,setm->Nstate,setm->Nstate);
                transpose_conjugate(sets->U_d2a,tempcm1,setm->Nstate);
                cc_matmul(tempcm1,sets->U_ref,overlap,setm->Nstate,setm->Nstate,setm->Nstate);

                memset(overlap2, 0, setm->Nstate * setm->Nstate * sizeof(double));
            
                for (int i = 0; i < setm->Nstate * setm->Nstate; i++){
                    tempdm1[i] = cabs(overlap[i]);
                }

                for (int i = 0; i < setm->Nstate; i++) {
                    idloc=maxloc(tempdm1, setm->Nstate * setm->Nstate);
                    id1 = idloc / setm->Nstate; 
                    id2 = idloc % setm->Nstate; 
                    overlap2[idloc] = (creal(overlap[idloc]) >= 0.0) ? 1.0 : -1.0;
                    for (j = 0; j < setm->Nstate; j++) {
                        tempdm1[id1 * setm->Nstate + j] = 0;
                        tempdm1[j * setm->Nstate + id2] = 0;
                    }
                }

                // printf("ovERLAP:\n");
                // for (int i = 0; i < setm->Nstate; i++) {
                //     for (int j = 0; j < setm->Nstate; j++) {
                //         printf("%f  ", overlap2[i * setm->Nstate + j]);
                //     }
                //     printf("\n");
                // }

                cd_matmul(sets->U_d2a, overlap2, tempcm1, setm->Nstate, setm->Nstate, setm->Nstate);
                memcpy(sets->U_d2a,tempcm1,setm->Nstate * setm->Nstate * sizeof(double complex));

                for (int i = 0; i < setm->Nstate * setm->Nstate; i++){
                    tempdm2[i] = fabs(overlap2[i]);
                }
                dd_matmul(sets->E_adia, tempdm2, tempdv1, 1, setm->Nstate, setm->Nstate);
                memcpy(sets->E_adia, tempdv1, setm->Nstate * sizeof(double));

                memset(tempdm1, 0, setm->Nstate * setm->Nstate * sizeof(double));
                for (int i = 0; i < setm->Nstate; i++) {
                    tempdm1[i * setm->Nstate + i] = sets->E_adia[i];
                }
                transpose_conjugate(sets->U_d2a, tempcm1, setm->Nstate);
                dc_matmul(tempdm1, tempcm1, tempcm2, setm->Nstate, setm->Nstate, setm->Nstate);
                cc_matmul(sets->U_d2a,tempcm2,sets->V,setm->Nstate, setm->Nstate, setm->Nstate);

                // printf("V after corr:\n");
                // for (int i = 0; i < setm->Nstate; i++) {
                //     for (int j = 0; j < setm->Nstate; j++) {
                //         printf("%18.8E  ", creal(sets->V[i * setm->Nstate + j]));
                //     }
                //     printf("\n");
                // }
                // printf("==========================================\n");

            }
            memcpy(sets->U_d2a_old, sets->U_ref, setm->Nstate * setm->Nstate * sizeof(double complex));
            memcpy(sets->U_ref, sets->U_d2a, setm->Nstate * setm->Nstate * sizeof(double complex));
            sets->if_ad_nac = 1;
        }
        
    

}
#endif 
