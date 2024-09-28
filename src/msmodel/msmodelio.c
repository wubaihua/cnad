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


void readinp_SEM(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "k_eff_SEM")){
        list=cJSON_GetObjectItem(item, "k_eff_SEM");
        if (list->type == cJSON_Number) {
            setm->k_eff_SEM = list->valuedouble; 
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "bias_SEM")) {
        list = cJSON_GetObjectItem(item, "bias_SEM");
        if (list->type == cJSON_Number) {
            setm->bias_SEM = list->valuedouble;
            setm->bias_SEM /= au_2_eV;
        }
    }

    
    if (NULL != cJSON_GetObjectItem(item, "delta_SEM")) {
        list = cJSON_GetObjectItem(item, "delta_SEM");
        if (list->type == cJSON_Number) {
            setm->delta_SEM = list->valuedouble; 
            setm->delta_SEM /= au_2_eV;
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "omega_c_SEM")) {
        list = cJSON_GetObjectItem(item, "omega_c_SEM");
        if (list->type == cJSON_Number) {
            setm->omega_c_SEM = list->valuedouble; 
            setm->omega_c_SEM *= au_2_ps*10.0/53.09;
        }
    }

    if (NULL !=  cJSON_GetObjectItem(item, "N_bath_SEM")){
        list=cJSON_GetObjectItem(item, "N_bath_SEM");
        setm->N_bath_SEM = list->valueint; 
    }

    setm->lambda_SEM = setm->k_eff_SEM * setm->k_eff_SEM * setm->omega_c_SEM * 0.5 ;
    setm->Nstate_SEM = 2;

    *Ndof1 = setm->Nstate_SEM;
    *Ndof2 = setm->N_bath_SEM;
    *Nstate = setm->Nstate_SEM;



    setm->H_ele_SEM = (double *)malloc(setm->Nstate_SEM * setm->Nstate_SEM * sizeof(double));

    setm->H_ele_SEM[0]=setm->bias_SEM * 0.5;
    setm->H_ele_SEM[1]=setm->delta_SEM;
    setm->H_ele_SEM[2]=setm->delta_SEM;
    setm->H_ele_SEM[3]=setm->bias_SEM * (-0.5);

    //debug
    // printf("N_bath_SEM =  %d\n",  setm->N_bath_SEM );
    // printf("k_eff_SEM  =  %f\n",  setm->k_eff_SEM  );
    // printf("bias_SEM   =  %f\n",  setm->bias_SEM   );
    // printf("delta_SEM  =  %f\n",  setm->delta_SEM  );
    // printf("lambda_SEM =  %f\n",  setm->lambda_SEM );
    // printf("omega_c_SEM=  %f\n",  setm->omega_c_SEM);
    // printf("F: %d\n", setm->Nstate_SEM);

    //debug
}



void readinp_msmodel(cJSON *json, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    if (strcmp(setm->msmodelname, "SBM") == 0 ||
       strcmp(setm->msmodelname, "sbm") == 0) {
        readinp_SBM(json, Ndof1, Ndof2, Nstate, setm);
    } else if (strcmp(setm->msmodelname, "morse3") == 0 ||
       strcmp(setm->msmodelname, "Morse3") == 0) {
        readinp_morse3(json, Ndof1, Ndof2, Nstate, setm);
    } else if (strcmp(setm->msmodelname, "SEM") == 0 ||
       strcmp(setm->msmodelname, "sem") == 0) {
        readinp_SEM(json, Ndof1, Ndof2, Nstate, setm);
         
    }


}

