

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
void parameter_SBM(double *mass, struct set_host *setm) {

    for (int j = 0; j < setm->N_bath_SBM; j++) {
        mass[j] = 1.0;
    }
    setm->c_SBM = (double *)malloc(setm->N_bath_SBM*sizeof(double));
    setm->omega_SBM = (double *)malloc(setm->N_bath_SBM*sizeof(double));
    switch (setm->bathtype) {
        case 1: // Ohmic, para=alpha
            for (int j = 0; j < setm->N_bath_SBM; j++) {
                setm->omega_SBM[j] = -setm->omega_c_SBM * log(1 - (double)(j + 1) / (setm->N_bath_SBM + 1));
                setm->c_SBM[j] = setm->omega_SBM[j] * sqrt(setm->alpha_SBM * setm->omega_c_SBM / (setm->N_bath_SBM + 1));
            }
            break;
        case 2: // Debye, para=lambda
            for (int j = 0; j < setm->N_bath_SBM; j++) {
                setm->omega_SBM[j] = setm->omega_c_SBM * tan(0.5 * M_PI * (1 - (double)(j + 1) / (setm->N_bath_SBM + 1)));
                setm->c_SBM[j] = setm->omega_SBM[j] * sqrt(setm->lambda_SBM * 2 / (setm->N_bath_SBM + 1));
            }
            break;
        case 3: // sub-Ohmic, para=alpha, s
            for (int j = 0; j < setm->N_bath_SBM; j++) {
                setm->omega_SBM[j] = -setm->omega_c_SBM * log(1 - (double)(j + 1) / (setm->N_bath_SBM + 1));
                setm->c_SBM[j] = sqrt(setm->alpha_SBM * pow(setm->omega_SBM[j], 1.0 + setm->s_SBM) * pow(setm->omega_c_SBM, 2.0 - setm->s_SBM) / (setm->N_bath_SBM + 1));
            }
            break;
        default:
            fprintf(stderr, "ERROR: unknown Spin-Boson Model bath kind\n");
            exit(1);
    }
}

// Sample the initial conditionals for trajectories of the model
void sample_SBM(double *P, double *R, double beta, struct set_host *setm) {
    
    double x2;
    for (int j = 0; j < setm->N_bath_SBM; j++) {
        if (beta > 99999) {
            if(setm->if_classical == 0){
                box_muller(&P[j], &x2, sqrt(0.5 * hbar * setm->omega_SBM[j]), 0.0);
                box_muller(&R[j], &x2, sqrt(0.5 * hbar / setm->omega_SBM[j]), 0.0);
            } else {
                P[j] = 0;
                R[j] = 0;
            }
            
        } else {
            if(setm->if_classical == 0){
                box_muller(&P[j], &x2, sqrt(0.5 * hbar * setm->omega_SBM[j] / tanh(0.5 * beta * hbar * setm->omega_SBM[j])), 0.0);
                box_muller(&R[j], &x2, sqrt(0.5 * hbar / (tanh(0.5 * beta * hbar * setm->omega_SBM[j]) * setm->omega_SBM[j])), 0.0);
            // P[j]=1,R[j]=1;   //debug 
            } else {
                box_muller(&P[j], &x2, sqrt(1.0 / beta), 0.0);
                box_muller(&R[j], &x2, sqrt(1.0 / (beta * setm->omega_SBM[j] * setm->omega_SBM[j])), 0.0);
            }
        
        }
    }
}

// Build the diabatic potential matrix of the model
void V_SBM(double *R, double *H, int forcetype, struct set_host *setm) {
    double sum1=0, sum2=0;
    for(int i=0;i<setm->N_bath_SBM;i++){
        sum1 += setm->c_SBM[i] * R[i];
        sum2 += setm->omega_SBM[i] * setm->omega_SBM[i] * R[i] * R[i];
    }
    switch (forcetype) {
        case 0:
            H[0] = setm->eps_SBM + sum1 + 0.5 * sum2;
            H[1] = setm->delta_SBM;
            H[3] = -setm->eps_SBM - sum1 + 0.5 * sum2;
            H[2] = setm->delta_SBM;
            break;
        case 1:
            H[0] = setm->eps_SBM + sum1;
            H[1] = setm->delta_SBM;
            H[3] = -setm->eps_SBM - sum1;
            H[2] = setm->delta_SBM;
            break;
        default:
            fprintf(stderr, "ERROR: unknown forcetype\n");
            exit(1);
    }

    // int slavecore_id=athread_get_id(-1);
//     if(slavecore_id == 0){
//         printf("c=%18.8E %18.8E w=%18.8E %18.8E\n", setm->c_SBM[0],setm->c_SBM[1],setm->omega_SBM[0],setm->omega_SBM[1]);
//         printf("R=%18.8E %18.8E V=%18.8E %18.8E %18.8E %18.8E\n", R[0],R[1],H[0],H[1],H[2],H[3]);
       
//     }
}

// Build the first-order derivative matrix of the model
void dV_SBM(double *R, double *dH, int forcetype, struct set_host *setm) {
    // for (int i = 0; i < 2; i++) {
    //     for (int j = 0; j < 2; j++) {
    //         for (int k = 0; k < N_bath_SBM; k++) {
    //             dH[i*2*N_bath_SBM+j*N_bath_SBM+k] = 0.0;
    //         }
    //     }
    // }
    memset(dH,0,4*setm->N_bath_SBM*sizeof(double));

    switch (forcetype) {
        case 0:
            for (int j = 0; j < setm->N_bath_SBM; j++) {
                dH[0*2*setm->N_bath_SBM+0*setm->N_bath_SBM+j] = setm->c_SBM[j]+setm->omega_SBM[j]* setm->omega_SBM[j] * R[j];
                dH[1*2*setm->N_bath_SBM+1*setm->N_bath_SBM+j] = -setm->c_SBM[j]+setm->omega_SBM[j]* setm->omega_SBM[j] * R[j];
              
            }
            break;
        case 1:
            for (int j = 0; j < setm->N_bath_SBM; j++) {
                dH[0*2*setm->N_bath_SBM+0*setm->N_bath_SBM+j] = setm->c_SBM[j];
                dH[1*2*setm->N_bath_SBM+1*setm->N_bath_SBM+j] = -setm->c_SBM[j];
            }
            break;
        default:
            fprintf(stderr, "ERROR: unknown forcetype\n");
            exit(1);
    }
}

// Calculate the nuclear force of the model
void nucforce_SBM(double *R, double *nf, struct set_host *setm) {
    for (int j = 0; j < setm->N_bath_SBM; j++) {
        nf[j] = setm->omega_SBM[j]* setm->omega_SBM[j] * R[j];
    }
    
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