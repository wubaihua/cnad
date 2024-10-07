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



    // setm->H_ele_SEM = (double *)malloc(setm->Nstate_SEM * setm->Nstate_SEM * sizeof(double));

    // setm->H_ele_SEM[0]=setm->bias_SEM * 0.5;
    // setm->H_ele_SEM[1]=setm->delta_SEM;
    // setm->H_ele_SEM[2]=setm->delta_SEM;
    // setm->H_ele_SEM[3]=setm->bias_SEM * (-0.5);

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



void readinp_FMO(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "lambda_SEM")){
        list=cJSON_GetObjectItem(item, "lambda_SEM");
        if (list->type == cJSON_Number) {
            setm->lambda_SEM = list->valuedouble; 
            setm->lambda_SEM *= au_2_ps*10.0/53.09;
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
    setm->Nstate_SEM = 7;

    *Ndof1 = setm->Nstate_SEM;
    *Ndof2 = setm->N_bath_SEM;
    *Nstate = setm->Nstate_SEM;



    // setm->H_ele_SEM = (double *)malloc(setm->Nstate_SEM * setm->Nstate_SEM * sizeof(double));

    
}


void readinp_SF(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "lambda_SEM")){
        list=cJSON_GetObjectItem(item, "lambda_SEM");
        if (list->type == cJSON_Number) {
            setm->lambda_SEM = list->valuedouble; 
            setm->lambda_SEM /= au_2_eV;
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "omega_c_SEM")) {
        list = cJSON_GetObjectItem(item, "omega_c_SEM");
        if (list->type == cJSON_Number) {
            setm->omega_c_SEM = list->valuedouble; 
            setm->omega_c_SEM /= au_2_eV;
        }
    }

    if (NULL !=  cJSON_GetObjectItem(item, "N_bath_SEM")){
        list=cJSON_GetObjectItem(item, "N_bath_SEM");
        setm->N_bath_SEM = list->valueint; 
    }
    setm->Nstate_SEM = 3;

    *Ndof1 = setm->Nstate_SEM;
    *Ndof2 = setm->N_bath_SEM;
    *Nstate = setm->Nstate_SEM;



    // setm->H_ele_SEM = (double *)malloc(setm->Nstate_SEM * setm->Nstate_SEM * sizeof(double));

    // setm->H_ele_SEM = (double *)malloc(setm->Nstate_SEM * setm->Nstate_SEM * sizeof(double));

    
}

void readinp_AIC(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "Nstate_AIC")){
        list=cJSON_GetObjectItem(item, "Nstate_AIC");
        setm->Nstate_aic = list->valueint; 
        
    }

    if (NULL !=  cJSON_GetObjectItem(item, "N_mode_AIC")){
        list=cJSON_GetObjectItem(item, "N_mode_AIC");
        setm->N_mode_aic = list->valueint; 
    }
   

    *Ndof1 = 1;
    *Ndof2 = setm->N_mode_aic;
    *Nstate = setm->Nstate_aic;
    
}



void readinp_pyrazine(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "type_pyrazine")){
        list=cJSON_GetObjectItem(item, "type_pyrazine");
        setm->type_pyrazine = list->valueint;         
    }

    switch (setm->type_pyrazine) {
        case 1:
            setm->Nstate_lvcm = 2;
            setm->N_mode_lvcm = 3;
            break;
        case 2:
            setm->Nstate_lvcm = 2;
            setm->N_mode_lvcm = 24;
            break;
    }

    *Ndof1 = 1;
    *Ndof2 = setm->N_mode_lvcm;
    *Nstate = setm->Nstate_lvcm;
    
}


void readinp_crco5(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "type_crco5")){
        list=cJSON_GetObjectItem(item, "type_crco5");
        setm->type_crco5 = list->valueint;         
    }

    switch (setm->type_crco5) {
        case 1:
            setm->Nstate_lvcm = 3;
            setm->N_mode_lvcm = 2;
            break;
        case 2:
            setm->Nstate_lvcm = 3;
            setm->N_mode_lvcm = 5;
            break;
    }

    *Ndof1 = 1;
    *Ndof2 = setm->N_mode_lvcm;
    *Nstate = setm->Nstate_lvcm;

}

void readinp_tully(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "type_tully")){
        list=cJSON_GetObjectItem(item, "type_tully");
        setm->type_tully = list->valueint;         
    }

    if (NULL !=  cJSON_GetObjectItem(item, "P0_tully")){
        list=cJSON_GetObjectItem(item, "P0_tully");
        if (list->type == cJSON_Number) {
            setm->P0_tully = list->valuedouble; 
        }        
    }

    
    *Ndof1 = 1;
    *Ndof2 = 1;
    *Nstate = 2;
    
}

void readinp_SEMdp(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    

    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "type_cal_SEMdp")){
        list=cJSON_GetObjectItem(item, "type_cal_SEMdp");
        
        setm->type_cal_SEMdp = list->valueint; 
        
    }

    if (NULL != cJSON_GetObjectItem(item, "bias_SEMdp")) {
        list = cJSON_GetObjectItem(item, "bias_SEMdp");
        if (list->type == cJSON_Number) {
            setm->bias_SEMdp = list->valuedouble;
            setm->bias_SEMdp /= au_2_wn;
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "lambda_SEMdp")) {
        list = cJSON_GetObjectItem(item, "lambda_SEMdp");
        if (list->type == cJSON_Number) {
            setm->lambda_SEMdp = list->valuedouble; 
            setm->delta_SEMdp /= au_2_wn;
        }
    }

    
    if (NULL != cJSON_GetObjectItem(item, "delta_SEMdp")) {
        list = cJSON_GetObjectItem(item, "delta_SEMdp");
        if (list->type == cJSON_Number) {
            setm->delta_SEMdp = list->valuedouble; 
            setm->delta_SEMdp /= au_2_wn;
        }
    }

    if (NULL != cJSON_GetObjectItem(item, "omega_c_SEMdp")) {
        list = cJSON_GetObjectItem(item, "omega_c_SEMdp");
        if (list->type == cJSON_Number) {
            setm->omega_c_SEMdp = list->valuedouble; 
            setm->omega_c_SEMdp *= au_2_ps*10.0/53.09;
        }
    }

    if (NULL !=  cJSON_GetObjectItem(item, "N_bath_SEMdp")){
        list=cJSON_GetObjectItem(item, "N_bath_SEMdp");
        setm->N_bath_SEMdp = list->valueint; 
    }

   
    setm->Nstate_SEMdp = 3;

    *Ndof1 = setm->Nstate_SEMdp - 1;
    *Ndof2 = setm->N_bath_SEMdp;
    *Nstate = setm->Nstate_SEMdp;



    // setm->H_ele_SEMdp = (double *)malloc(setm->Nstate_SEMdp * setm->Nstate_SEMdp * sizeof(double));

    // for (int i = 0; i < 3 * 3; i++) setm->H_ele_SEMdp[i] = 0.0;
    // setm->H_ele_SEMdp[0 * 3 + 0] = setm->bias_SEMdp + 1000.0 / au_2_wn;
    // setm->H_ele_SEMdp[0 * 3 + 1] = setm->delta_SEMdp;
    // setm->H_ele_SEMdp[1 * 3 + 1] = 1000.0 / au_2_wn;
    // setm->H_ele_SEMdp[1 * 3 + 0] = setm->delta_SEMdp;

    // if (setm->type_cal_SEMdp == 2) {
    //     setm->H_ele_SEMdp[0 * 3 + 2] = 2000 * au_2_ps * 10.0 / 53.09;
    //     setm->H_ele_SEMdp[2 * 3 + 0] = 2000 * au_2_ps * 10.0 / 53.09;
    //     setm->H_ele_SEMdp[1 * 3 + 2] = -400 * au_2_ps * 10.0 / 53.09;
    //     setm->H_ele_SEMdp[2 * 3 + 1] = -400 * au_2_ps * 10.0 / 53.09;
    // }

    // setm->dipole_SEMdp = (double *)malloc(3 * sizeof(double));
    // setm->dipole_SEMdp[0] = 5.0;
    // setm->dipole_SEMdp[1] = -1.0;
    // setm->dipole_SEMdp[2] = 0.0;

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


void readinp_FMOdp(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate, struct set_host *setm) {
    cJSON *list;
    
    if (NULL !=  cJSON_GetObjectItem(item, "type_cal_SEMdp")){
        list=cJSON_GetObjectItem(item, "type_cal_SEMdp");
        
        setm->type_cal_SEMdp = list->valueint; 
        
    }

    if (NULL != cJSON_GetObjectItem(item, "lambda_SEMdp")) {
        list = cJSON_GetObjectItem(item, "lambda_SEMdp");
        if (list->type == cJSON_Number) {
            setm->lambda_SEMdp = list->valuedouble; 
            setm->delta_SEMdp /= au_2_wn;
        }
    }

  

    if (NULL != cJSON_GetObjectItem(item, "omega_c_SEMdp")) {
        list = cJSON_GetObjectItem(item, "omega_c_SEMdp");
        if (list->type == cJSON_Number) {
            setm->omega_c_SEMdp = list->valuedouble; 
            setm->omega_c_SEMdp *= au_2_ps*10.0/53.09;
        }
    }

    if (NULL !=  cJSON_GetObjectItem(item, "N_bath_SEMdp")){
        list=cJSON_GetObjectItem(item, "N_bath_SEMdp");
        setm->N_bath_SEMdp = list->valueint; 
    }
    
    setm->Nstate_SEMdp = 8;
    *Nstate = setm->Nstate_SEMdp;
    *Ndof1 = setm->Nstate_SEMdp - 1;
    *Ndof2 = setm->N_bath_SEMdp;

    // setm->H_ele_SEMdp = (double *)malloc(8 * 8 * sizeof(double));
    // for (int i = 0; i < 8 * 8; i++) setm->H_ele_SEMdp[i] = 0.0;

    // double H_ele_values[7][7] = {
    //     {12410.0, -87.7, 5.5, -5.9, 6.7, -13.7, -9.9},
    //     {-87.7, 12530.0, 30.8, 8.2, 0.7, 11.8, 4.3},
    //     {5.5, 30.8, 12210.0, -53.5, -2.2, -9.6, 6.0},
    //     {-5.9, 8.2, -53.5, 12320.0, -70.7, -17.0, -63.3},
    //     {6.7, 0.7, -2.2, -70.7, 12480.0, 81.1, -1.3},
    //     {-13.7, 11.8, -9.6, -17.0, 81.1, 12630.0, 39.7},
    //     {-9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 12440.0}
    // };

    // for (int i = 0; i < 7; i++) {
    //     for (int j = 0; j < 7; j++) {
    //         setm->H_ele_SEMdp[i * 8 + j] = H_ele_values[i][j] / au_2_wn;
    //     }
    // }

   
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
         
    } else if (strcmp(setm->msmodelname, "FMO") == 0 ||
       strcmp(setm->msmodelname, "fmo") == 0) {
        readinp_FMO(json, Ndof1, Ndof2, Nstate, setm);
         
    } else if (strcmp(setm->msmodelname, "SF") == 0 ||
       strcmp(setm->msmodelname, "sf") == 0) {
        readinp_SF(json, Ndof1, Ndof2, Nstate, setm);
         
    } else if (strcmp(setm->msmodelname, "AIC") == 0 ||
       strcmp(setm->msmodelname, "aic") == 0) {
        readinp_AIC(json, Ndof1, Ndof2, Nstate, setm);
         
    } else if (strcmp(setm->msmodelname, "pyrazine") == 0) {
        readinp_pyrazine(json, Ndof1, Ndof2, Nstate, setm);
         
    } else if (strcmp(setm->msmodelname, "crco5") == 0 ||
       strcmp(setm->msmodelname, "CrCO5") == 0) {
        readinp_crco5(json, Ndof1, Ndof2, Nstate, setm);
         
    } else if (strcmp(setm->msmodelname, "tully") == 0) {
        readinp_tully(json, Ndof1, Ndof2, Nstate, setm);

    } else if (strcmp(setm->msmodelname, "SEMdp") == 0 ||
        strcmp(setm->msmodelname, "semdp") == 0) {
        readinp_SEMdp(json, Ndof1, Ndof2, Nstate, setm);

    } else if (strcmp(setm->msmodelname, "FMOdp") == 0 ||
        strcmp(setm->msmodelname, "fmodp") == 0) {
        readinp_FMOdp(json, Ndof1, Ndof2, Nstate, setm);
    }



}

