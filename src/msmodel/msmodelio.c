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
#include "msmodelio.h"
#include "def_host.h"



// int forcetype;
// char msmodelname[200];


// // Spin-Boson Model parameters
// int N_bath_SBM, bathtype; // bathtype=1 for Ohmic; bathtype=2 for Debye
// double eps_SBM, delta_SBM, alpha_SBM, omega_c_SBM, lambda_SBM, s_SBM;


void readinp_SBM(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "N_bath_SBM")){
        list=cJSON_GetObjectItem(item, "N_bath_SBM");
        setm->N_bath_SBM = list->valueint; 
    }

    if (NULL != cJSON_GetObjectItem(item, "bathtype")) {
        list = cJSON_GetObjectItem(item, "bathtype");
        setm->bathtype = list->valueint; 
    }

    if (NULL != cJSON_GetObjectItem(item, "eps_SBM")) {
        list = cJSON_GetObjectItem(item, "eps_SBM");
        if (list->type == cJSON_Number) {
            setm->eps_SBM = list->valuedouble;
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "delta_SBM")) {
        list = cJSON_GetObjectItem(item, "delta_SBM");
        if (list->type == cJSON_Number) {
            setm->delta_SBM = list->valuedouble; 
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "alpha_SBM")) {
        list = cJSON_GetObjectItem(item, "alpha_SBM");
        if (list->type == cJSON_Number) {
            setm->alpha_SBM = list->valuedouble;
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "omega_c_SBM")) {
        list = cJSON_GetObjectItem(item, "omega_c_SBM");
        if (list->type == cJSON_Number) {
            setm->omega_c_SBM = list->valuedouble; 
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "lambda_SBM")) {
        list = cJSON_GetObjectItem(item, "lambda_SBM");
        if (list->type == cJSON_Number) {
            setm->lambda_SBM = list->valuedouble; 
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "s_SBM")) {
        list = cJSON_GetObjectItem(item, "s_SBM");
        if (list->type == cJSON_Number) {
            setm->s_SBM = list->valuedouble; 
        }
    }

    


    *Ndof1 = 1;
    *Ndof2 = setm->N_bath_SBM;
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



void readinp_morse3(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "type_morse3")){
        list=cJSON_GetObjectItem(item, "type_morse3");
        setm->type_morse3 = list->valueint; 
    }

    


    *Ndof1 = 1;
    *Ndof2 = 1;
    *Nstate = 3;


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






void readinp_msmodel(cJSON *json, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        readinp_SBM(json, Ndof1, Ndof2, Nstate, setm);
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        readinp_morse3(json, Ndof1, Ndof2, Nstate, setm);
    } 


}

