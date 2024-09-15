

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



// Spin-Boson Model parameters
// int N_bath_SBM, bathtype; // bathtype=1 for Ohmic; bathtype=2 for Debye
// double eps_SBM, delta_SBM, alpha_SBM, omega_c_SBM, lambda_SBM, s_SBM;
double *c_SBM, *omega_SBM;

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
void parameter_SBM(double *mass) {

    for (int j = 0; j < N_bath_SBM; j++) {
        mass[j] = 1.0;
    }
    c_SBM = (double *)malloc(N_bath_SBM*sizeof(double));
    omega_SBM = (double *)malloc(N_bath_SBM*sizeof(double));
    switch (bathtype) {
        case 1: // Ohmic, para=alpha
            for (int j = 0; j < N_bath_SBM; j++) {
                omega_SBM[j] = -omega_c_SBM * log(1 - (double)(j + 1) / (N_bath_SBM + 1));
                c_SBM[j] = omega_SBM[j] * sqrt(alpha_SBM * omega_c_SBM / (N_bath_SBM + 1));
            }
            break;
        case 2: // Debye, para=lambda
            for (int j = 0; j < N_bath_SBM; j++) {
                omega_SBM[j] = omega_c_SBM * tan(0.5 * pi * (1 - (double)(j + 1) / (N_bath_SBM + 1)));
                c_SBM[j] = omega_SBM[j] * sqrt(lambda_SBM * 2 / (N_bath_SBM + 1));
            }
            break;
        case 3: // sub-Ohmic, para=alpha, s
            for (int j = 0; j < N_bath_SBM; j++) {
                omega_SBM[j] = -omega_c_SBM * log(1 - (double)(j + 1) / (N_bath_SBM + 1));
                c_SBM[j] = sqrt(alpha_SBM * pow(omega_SBM[j], 1.0 + s_SBM) * pow(omega_c_SBM, 2.0 - s_SBM) / (N_bath_SBM + 1));
            }
            break;
        default:
            fprintf(stderr, "ERROR: unknown Spin-Boson Model bath kind\n");
            exit(1);
    }
}

// Sample the initial conditionals for trajectories of the model
void sample_SBM(double *P, double *R, double beta) {
    
    double x2;
    for (int j = 0; j < N_bath_SBM; j++) {
        if (beta > 99999) {
            box_muller(&P[j], &x2, sqrt(0.5 * hbar * omega_SBM[j]), 0.0);
            box_muller(&R[j], &x2, sqrt(0.5 * hbar / omega_SBM[j]), 0.0);
        } else {
            box_muller(&P[j], &x2, sqrt(0.5 * hbar * omega_SBM[j] / tanh(0.5 * beta * hbar * omega_SBM[j])), 0.0);
            box_muller(&R[j], &x2, sqrt(0.5 * hbar / (tanh(0.5 * beta * hbar * omega_SBM[j]) * omega_SBM[j])), 0.0);
            // P[j]=0,R[j]=0;   //debug 
        
        }
    }
}

// Build the diabatic potential matrix of the model
void V_SBM(double *R, double *H, int forcetype) {
    double sum1=0, sum2=0;
    for(int i=0;i<N_bath_SBM;i++){
        sum1 += c_SBM[i] * R[i];
        sum2 += omega_SBM[i] * omega_SBM[i] + R[i] * R[i];
    }
    switch (forcetype) {
        case 0:
            H[0] = eps_SBM + sum1 + 0.5 * sum2;
            H[1] = delta_SBM;
            H[3] = -eps_SBM - sum1 + 0.5 * sum2;
            H[2] = delta_SBM;
            break;
        case 1:
            H[0] = eps_SBM + sum1;
            H[1] = delta_SBM;
            H[3] = -eps_SBM - sum1;
            H[2] = delta_SBM;
            break;
        default:
            fprintf(stderr, "ERROR: unknown forcetype\n");
            exit(1);
    }
}

// Build the first-order derivative matrix of the model
void dV_SBM(double *R, double *dH, int forcetype) {
    // for (int i = 0; i < 2; i++) {
    //     for (int j = 0; j < 2; j++) {
    //         for (int k = 0; k < N_bath_SBM; k++) {
    //             dH[i*2*N_bath_SBM+j*N_bath_SBM+k] = 0.0;
    //         }
    //     }
    // }
    memset(dH,0,4*N_bath_SBM*sizeof(double));

    switch (forcetype) {
        case 0:
            for (int j = 0; j < N_bath_SBM; j++) {
                dH[0*2*N_bath_SBM+0*N_bath_SBM+j] = c_SBM[j]+omega_SBM[j]* omega_SBM[j] * R[j];
                dH[1*2*N_bath_SBM+1*N_bath_SBM+j] = -c_SBM[j]+omega_SBM[j]* omega_SBM[j] * R[j];
              
            }
            break;
        case 1:
            for (int j = 0; j < N_bath_SBM; j++) {
                dH[0*2*N_bath_SBM+0*N_bath_SBM+j] = c_SBM[j];
                dH[1*2*N_bath_SBM+1*N_bath_SBM+j] = -c_SBM[j];
            }
            break;
        default:
            fprintf(stderr, "ERROR: unknown forcetype\n");
            exit(1);
    }
}

// Calculate the nuclear force of the model
void nucforce_SBM(double *R, double *nf) {
    for (int j = 0; j < N_bath_SBM; j++) {
        nf[j] = omega_SBM[j]* omega_SBM[j] * R[j];
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