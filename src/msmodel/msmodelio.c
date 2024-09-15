#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"  // Assuming you have a math.h header for functions like box_muller
// #include "constant.h"
#include <string.h>
#include <stdio.h>
#include "cJSON.h"


int forcetype;
char msmodelname[200];


// Spin-Boson Model parameters
int N_bath_SBM, bathtype; // bathtype=1 for Ohmic; bathtype=2 for Debye
double eps_SBM, delta_SBM, alpha_SBM, omega_c_SBM, lambda_SBM, s_SBM;
void readinp_SBM(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "N_bath_SBM")){
            list=cJSON_GetObjectItem(item, "N_bath_SBM");
            N_bath_SBM = list->valueint; 
    }

    if (NULL != cJSON_GetObjectItem(item, "bathtype")) {
    list = cJSON_GetObjectItem(item, "bathtype");
    bathtype = list->valueint; 
    }

    if (NULL != cJSON_GetObjectItem(item, "eps_SBM")) {
        list = cJSON_GetObjectItem(item, "eps_SBM");
        if (list->type == cJSON_Number) {
            eps_SBM = list->valuedouble;
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "delta_SBM")) {
        list = cJSON_GetObjectItem(item, "delta_SBM");
        if (list->type == cJSON_Number) {
            delta_SBM = list->valuedouble; 
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "alpha_SBM")) {
        list = cJSON_GetObjectItem(item, "alpha_SBM");
        if (list->type == cJSON_Number) {
            alpha_SBM = list->valuedouble;
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "omega_c_SBM")) {
        list = cJSON_GetObjectItem(item, "omega_c_SBM");
        if (list->type == cJSON_Number) {
            omega_c_SBM = list->valuedouble; 
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "lambda_SBM")) {
        list = cJSON_GetObjectItem(item, "lambda_SBM");
        if (list->type == cJSON_Number) {
            lambda_SBM = list->valuedouble; 
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "s_SBM")) {
        list = cJSON_GetObjectItem(item, "s_SBM");
        if (list->type == cJSON_Number) {
            s_SBM = list->valuedouble; 
        }
    }

    


    *Ndof1 = 1;
    *Ndof2 = N_bath_SBM;
    *Nstate = 2;


    // //debug
    // printf("N_bath_SBM: %d\n", N_bath_SBM);
    // printf("bathtype: %d\n", bathtype);
    // printf("eps_SBM: %f\n", eps_SBM);
    // printf("delta_SBM: %f\n", delta_SBM);
    // printf("alpha_SBM: %f\n", alpha_SBM);
    // printf("omega_c_SBM: %f\n", omega_c_SBM);
    // printf("F: %d\n", Nstate);

    // //debug
}








void readinp_msmodel(cJSON *json, int *Ndof1, int *Ndof2, int *Nstate) {
    if (strcmp(trim(adjustl(msmodelname)), "SBM") == 0 ||
       strcmp(trim(adjustl(msmodelname)), "sbm") == 0) {
        readinp_SBM(json, Ndof1, Ndof2, Nstate);
    }
}

