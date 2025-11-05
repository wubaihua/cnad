

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
#include "sbm.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif


// Spin-Boson Model parameters
// int N_bath_SBM, bathtype; // bathtype=1 for Ohmic; bathtype=2 for Debye
// double eps_SBM, delta_SBM, alpha_SBM, omega_c_SBM, lambda_SBM, s_SBM;
// double *c_SBM, *omega_SBM;

// Read the model type from the input file
// void readinp_SBM(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate) {
    

//     cJSON *list;
    
//     if (NULL !=  cJSON_GetObjectItem(item, "N_bath_SBM")){
//             list=cJSON_GetObjectItem(item, "N_bath_SBM");
//             N_bath_SBM = list->valueint; 
//     }

//     if (NULL != cJSON_GetObjectItem(item, "bathtype")) {
//     list = cJSON_GetObjectItem(item, "bathtype");
//     bathtype = list->valueint; 
//     }

//     if (NULL != cJSON_GetObjectItem(item, "eps_SBM")) {
//         list = cJSON_GetObjectItem(item, "eps_SBM");
//         if (list->type == cJSON_Number) {
//             eps_SBM = list->valuedouble;
//         }
//     }

//     if (NULL != cJSON_GetObjectItem(item, "delta_SBM")) {
//         list = cJSON_GetObjectItem(item, "delta_SBM");
//         if (list->type == cJSON_Number) {
//             delta_SBM = list->valuedouble; 
//         }
//     }

//     if (NULL != cJSON_GetObjectItem(item, "alpha_SBM")) {
//         list = cJSON_GetObjectItem(item, "alpha_SBM");
//         if (list->type == cJSON_Number) {
//             alpha_SBM = list->valuedouble;
//         }
//     }

//     if (NULL != cJSON_GetObjectItem(item, "omega_c_SBM")) {
//         list = cJSON_GetObjectItem(item, "omega_c_SBM");
//         if (list->type == cJSON_Number) {
//             omega_c_SBM = list->valuedouble; 
//         }
//     }

    


//     *Ndof1 = 1;
//     *Ndof2 = N_bath_SBM;
//     *Nstate = 2;


//     // //debug
//     // printf("N_bath_SBM: %d\n", N_bath_SBM);
//     // printf("bathtype: %d\n", bathtype);
//     // printf("eps_SBM: %f\n", eps_SBM);
//     // printf("delta_SBM: %f\n", delta_SBM);
//     // printf("alpha_SBM: %f\n", alpha_SBM);
//     // printf("omega_c_SBM: %f\n", omega_c_SBM);
//     // printf("F: %d\n", Nstate);

//     // //debug
// }

// Initialize model parameters
void parameter_SBMTD(double *mass, struct set_host *setm) {
    setm->eps_SBM = 0.0;
    parameter_SBM(mass, setm);
    
}

// Sample the initial conditionals for trajectories of the model
void sample_SBMTD(double *P, double *R, double beta, struct set_host *setm) {
    
    sample_SBM(P, R, beta, setm);
}

// Build the diabatic potential matrix of the model
void V_SBMTD(double *R, double *H, int forcetype, double t, struct set_host *setm) {
    V_SBM(R, H, forcetype, setm);
    H[0] += setm->eps0_SBMTD * (setm->alpha_SBMTD + (1.0 - setm->alpha_SBMTD) * cos(setm->omega_SBMTD * t));
    H[3] -= setm->eps0_SBMTD * (setm->alpha_SBMTD + (1.0 - setm->alpha_SBMTD) * cos(setm->omega_SBMTD * t));
    
    
}

// Build the first-order derivative matrix of the model
void dV_SBMTD(double *R, double *dH, int forcetype, double t, struct set_host *setm) {
    dV_SBM(R, dH, forcetype, setm);

}

// Calculate the nuclear force of the model
void nucforce_SBMTD(double *R, double *nf, struct set_host *setm) {
    nucforce_SBM(R, nf, setm);
    
}

// Compute the cfweight of the model
// void cfweight_SBM(double w0[2][2], double wt[2][2], double beta) {
//     double rho[2][2], Heff[2][2];
//     double E[2], C[2][2], expe[2][2], F[2][2];
//     double lbd_ohmic;
    
//     for (int i = 0; i < 2; i++) {
//         for (int j = 0; j < 2; j++) {
//             Heff[i][j] = 0.0;
//             rho[i][j] = 0.0;
//             expe[i][j] = 0.0;
//             F[i][j] = 0.0;
//         }
//     }

//     Heff[0][1] = delta_SBM;
//     Heff[1][0] = delta_SBM;

//     lbd_ohmic = 0.5 * alpha_SBM * omega_c_SBM / 2;
//     for (int i = 0; i < 2; i++) {
//         for (int j = 0; j < 2; j++) {
//             Heff[i][j] = Heff[i][j] * exp(-lbd_ohmic * beta / 3);
//         }
//     }
// }