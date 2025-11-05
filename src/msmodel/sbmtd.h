

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
void parameter_SBMTD(double *mass, struct set_host *setm);
// Sample the initial conditionals for trajectories of the model
void sample_SBMTD(double *P, double *R, double beta, struct set_host *setm);
// Build the diabatic potential matrix of the model
void V_SBMTD(double *R, double *H, int forcetype, double t, struct set_host *setm);

// Build the first-order derivative matrix of the model
void dV_SBMTD(double *R, double *dH, int forcetype, double t, struct set_host *setm);

// Calculate the nuclear force of the model
void nucforce_SBMTD(double *R, double *nf, struct set_host *setm) ;

