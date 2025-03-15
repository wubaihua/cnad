#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"
#include <string.h>
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif
#include "lvcm.h"

// #define PI 3.141592653589793
// #define HBAR 1.0545718e-34

// int setm->N_mode_lvcm, setm->Nstate_lvcm;
// double setm->L_lvcm;
// double *setm->eps_lvcm, *setm->miu_lvcm, *setm->omega_lvcm, *setm->lambda_lvcm;

// void readinp_lvcm(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     fscanf(idinp, "%d", &setm->Nstate_lvcm);
//     fscanf(idinp, "%d", &setm->N_mode_lvcm);

//     *Ndof1 = 1;
//     *Ndof2 = setm->N_mode_lvcm;
//     *Nstate = setm->Nstate_lvcm;
// }


void parameter_crco5(double *mass, struct set_host *setm) {
    
    parameter_LVCM(mass, setm);



    switch (setm->type_crco5) {
        case 1:

            setm->H_ele_lvcm[0] = 0.0424  / au_2_eV ;
            setm->H_ele_lvcm[4] = 0.0424  / au_2_eV ;
            setm->H_ele_lvcm[8] = 0.4344  / au_2_eV ;
            setm->omega_lvcm[0] = 0.0129 / au_2_eV ;
            setm->omega_lvcm[1] = 0.0129 / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.0328  / au_2_eV ;
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0328 / au_2_eV ; 
    
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0328 / au_2_eV ; 
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0328 / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = -0.0978 / au_2_eV ; 
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = -0.0978 / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = -0.0978 / au_2_eV ; 
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.0978 / au_2_eV ;
            
            double R0[5]={0.0,14.3514,-9.9699,-7.0189,0.0};
            double alpha[5]={0.4501,0.4286,0.6204,0.4535,0.5539};
            for(int i = 0; i < setm->N_mode_lvcm; i++){
                setm->R0_lvcm[i] = R0[i] / sqrt(setm->omega_lvcm[i]);
                setm->alpha_lvcm[i] = alpha[i] * sqrt(2.0);
            }

            break;

        case 2:
            // setm->H_ele_lvcm[0] = 0.0424  / au_2_eV ;
            // setm->H_ele_lvcm[4] = 0.0424  / au_2_eV ;
            // setm->H_ele_lvcm[8] = 0.4344  / au_2_eV ;
            // setm->omega_lvcm[0] = 0.0129 / au_2_eV ;
            // setm->omega_lvcm[1] = 0.0129 / au_2_eV ;
            // setm->omega_lvcm[2] = 0.0342 / au_2_eV ;
            // setm->omega_lvcm[3] = 0.0561 / au_2_eV ;
            // setm->omega_lvcm[4] = 0.0561 / au_2_eV ;
            
            // setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.0328  / au_2_eV ;
            // setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0328 / au_2_eV ; 
    
            // setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0328 / au_2_eV ; 
            // setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0328 / au_2_eV ;
            
            // setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = -0.0978 / au_2_eV ; 
            // setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = -0.0978 / au_2_eV ;

            // setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = -0.0978 / au_2_eV ; 
            // setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.0978 / au_2_eV ;
            // double R0[5]={0.0,14.3514,-9.9699,-7.0189,0.0};
            // double alpha[5]={0.4501,0.4286,0.6204,0.4535,0.5539}
            // for(int i = 0; i < sethm->N_mode_lvcm; i++){
            //     setm->R0_lvcm[i] = R0[i] / sqrt(setm->omega_lvcm[i]);
            //     setm->alpha_lvcm[i] = alpha[i] * sqrt(2.0);
            // }

            break;
            
     }

     
}

void sample_crco5(double *P, double *R, struct set_host *setm) {
    int j;
    double x2;

    sample_LVCM(P, R, setm);
    

}

void V_crco5(double *R, double complex *H, int forcetype, struct set_host *setm) {
    
    V_LVCM(R, H, forcetype, setm);
    
}

void dV_crco5(double *R, double complex *dH, int forcetype, struct set_host *setm) {
    
    dV_LVCM(R, dH, forcetype, setm);

}

void nucforce_crco5(double *R, double *nf, struct set_host *setm) {

    nucforce_LVCM(R, nf, setm);

}



