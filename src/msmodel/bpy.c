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


void parameter_bpy(double *mass, struct set_host *setm) {
    
    parameter_LVCM(mass, setm);

    

    
    //ele: 0-2: T1, 3-5: T2, 6: S1, 7: S2, 8-10: T3
    //nuc: 0-3: 7(a'), 11(a'), 13(a'), 30(a'); 4-5: 8(a"), 23(a")
    setm->omega_lvcm[0] = 0.0116 / au_2_eV ;
    setm->omega_lvcm[1] = 0.0188 / au_2_eV ;
    setm->omega_lvcm[2] = 0.0229 / au_2_eV ;
    setm->omega_lvcm[3] = 0.0792 / au_2_eV ;
    setm->omega_lvcm[4] = 0.0118 / au_2_eV ;
    setm->omega_lvcm[5] = 0.0601 / au_2_eV ;  

    //W^{T1, T1}
    setm->H_ele_lvcm[0 * 11 + 0] = 2.81 / au_2_eV ;
    setm->H_ele_lvcm[1 * 11 + 1] = 2.81 / au_2_eV ;
    setm->H_ele_lvcm[2 * 11 + 2] = 2.81 / au_2_eV ;

    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.0161  / au_2_eV         ;
    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0002  / au_2_eV         ;
    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.0261 / au_2_eV         ;
    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = -0.0196 / au_2_eV         ;
        
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = -0.0161  / au_2_eV       ;
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0002  / au_2_eV       ;
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = -0.0261 / au_2_eV       ;
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = -0.0196 / au_2_eV       ;
        
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = -0.0161  / au_2_eV        ;
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0002  / au_2_eV        ;
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = -0.0261 / au_2_eV        ;
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = -0.0196 / au_2_eV        ;

    //W^{T2, T2}
    setm->H_ele_lvcm[3 * 11 + 3] = 2.93 / au_2_eV ;
    setm->H_ele_lvcm[4 * 11 + 4] = 2.93 / au_2_eV ;
    setm->H_ele_lvcm[5 * 11 + 5] = 2.93 / au_2_eV ;

    setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.019  / au_2_eV     ;
    setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = -0.0006  / au_2_eV   ;
    setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = -0.0322 / au_2_eV    ;
    setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0433 / au_2_eV     ;

    setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.019  / au_2_eV    ;
    setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = -0.0006  / au_2_eV  ;
    setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = -0.0322 / au_2_eV   ;
    setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0433 / au_2_eV    ;

    setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.019  / au_2_eV    ;
    setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = -0.0006  / au_2_eV  ;
    setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = -0.0322 / au_2_eV   ;
    setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0433 / au_2_eV    ;

    //W^{S1, S1}
    setm->H_ele_lvcm[6 * 11 + 6] = 2.94 / au_2_eV ;
    setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = -0.0172  / au_2_eV  ;
    setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.009  / au_2_eV     ;
    setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = -0.0289  / au_2_eV   ;
    setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = -0.0187  / au_2_eV    ;

    //W^{S2, S2}
    setm->H_ele_lvcm[7 * 11 + 7] = 3.11 / au_2_eV ;
    setm->c_lvcm[7 * setm->Nstate_lvcm * setm->N_mode_lvcm + 7 * setm->N_mode_lvcm + 0] = 0.0187  / au_2_eV    ;
    setm->c_lvcm[7 * setm->Nstate_lvcm * setm->N_mode_lvcm + 7 * setm->N_mode_lvcm + 1] = 0.0091  / au_2_eV    ;
    setm->c_lvcm[7 * setm->Nstate_lvcm * setm->N_mode_lvcm + 7 * setm->N_mode_lvcm + 2] = -0.0271  / au_2_eV   ;
    setm->c_lvcm[7 * setm->Nstate_lvcm * setm->N_mode_lvcm + 7 * setm->N_mode_lvcm + 3] = 0.0404  / au_2_eV    ;


    //W^{T3, T3}
    setm->H_ele_lvcm[8 * 11 + 8] =   3.22 / au_2_eV ;
    setm->H_ele_lvcm[9 * 11 + 9] = 3.22 / au_2_eV ;
    setm->H_ele_lvcm[10 * 11 + 10] = 3.22 / au_2_eV ;

    setm->c_lvcm[8 * setm->Nstate_lvcm * setm->N_mode_lvcm + 8 * setm->N_mode_lvcm + 0] = 0.0015  / au_2_eV     ;
    setm->c_lvcm[8 * setm->Nstate_lvcm * setm->N_mode_lvcm + 8 * setm->N_mode_lvcm + 1] = 0.0056  / au_2_eV     ;
    setm->c_lvcm[8 * setm->Nstate_lvcm * setm->N_mode_lvcm + 8 * setm->N_mode_lvcm + 2] = 0.0133 / au_2_eV      ;
    setm->c_lvcm[8 * setm->Nstate_lvcm * setm->N_mode_lvcm + 8 * setm->N_mode_lvcm + 3] = 0.0033 / au_2_eV      ;

    setm->c_lvcm[9 * setm->Nstate_lvcm * setm->N_mode_lvcm + 9 * setm->N_mode_lvcm + 0] = 0.0015  / au_2_eV     ;
    setm->c_lvcm[9 * setm->Nstate_lvcm * setm->N_mode_lvcm + 9 * setm->N_mode_lvcm + 1] = 0.0056  / au_2_eV     ;
    setm->c_lvcm[9 * setm->Nstate_lvcm * setm->N_mode_lvcm + 9 * setm->N_mode_lvcm + 2] = 0.0133 / au_2_eV     ;
    setm->c_lvcm[9 * setm->Nstate_lvcm * setm->N_mode_lvcm + 9 * setm->N_mode_lvcm + 3] = 0.0033 / au_2_eV     ;

    setm->c_lvcm[10 * setm->Nstate_lvcm * setm->N_mode_lvcm + 10 * setm->N_mode_lvcm + 0] = 0.0015  / au_2_eV    ;
    setm->c_lvcm[10 * setm->Nstate_lvcm * setm->N_mode_lvcm + 10 * setm->N_mode_lvcm + 1] = 0.0056  / au_2_eV    ;
    setm->c_lvcm[10 * setm->Nstate_lvcm * setm->N_mode_lvcm + 10 * setm->N_mode_lvcm + 2] = 0.0133 / au_2_eV    ;
    setm->c_lvcm[10 * setm->Nstate_lvcm * setm->N_mode_lvcm + 10 * setm->N_mode_lvcm + 3] = 0.0033 / au_2_eV    ;

    //W^{S1, S2}
    setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 7 * setm->N_mode_lvcm + 4] = 0.0114  / au_2_eV     ;
    setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 7 * setm->N_mode_lvcm + 5] = 0.0237  / au_2_eV     ;
    setm->c_lvcm[7 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0114  / au_2_eV     ;
    setm->c_lvcm[7 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0237  / au_2_eV     ;

    //W^{S2, T3}
    setm->H_ele_lvcm[7 * 11 + 8] = ( 0.0274 - I * 0.0056 ) / au_2_eV ;
    setm->H_ele_lvcm[7 * 11 + 10] = ( 0.0274 + I * 0.0056 ) / au_2_eV ;
    setm->H_ele_lvcm[8 * 11 + 7] = ( 0.0274 + I * 0.0056 ) / au_2_eV ;
    setm->H_ele_lvcm[10 * 11 + 7] = ( 0.0274 - I * 0.0056 ) / au_2_eV ;

    //W^{T1, S2}
    setm->H_ele_lvcm[0 * 11 + 7] = ( -0.0719 - I * 0.0196 ) / au_2_eV ;
    setm->H_ele_lvcm[2 * 11 + 7] = ( -0.0719 + I * 0.0196 ) / au_2_eV ;
    setm->H_ele_lvcm[7 * 11 + 0] = ( -0.0719 + I * 0.0196 ) / au_2_eV ;
    setm->H_ele_lvcm[7 * 11 + 2] = ( -0.0719 - I * 0.0196 ) / au_2_eV ;

    //W^{T2, S1}
    setm->H_ele_lvcm[3 * 11 + 6] = ( 0.0769 + I * 0.0186 ) / au_2_eV ;
    setm->H_ele_lvcm[5 * 11 + 6] = ( 0.0769 - I * 0.0186 ) / au_2_eV ;
    setm->H_ele_lvcm[6 * 11 + 3] = ( 0.0769 - I * 0.0186 ) / au_2_eV ;
    setm->H_ele_lvcm[6 * 11 + 5] = ( 0.0769 + I * 0.0186 ) / au_2_eV ;

    //W^{T1, T2}
    setm->H_ele_lvcm[1 * 11 + 3] = ( 0.0719 + I * 0.0177 ) / au_2_eV ;
    setm->H_ele_lvcm[0 * 11 + 4] = - ( 0.0719 - I * 0.0177 ) / au_2_eV ;
    setm->H_ele_lvcm[2 * 11 + 4] = ( 0.0719 + I * 0.0177 )  / au_2_eV ;
    setm->H_ele_lvcm[1 * 11 + 5] = - ( 0.0719 - I * 0.0177 ) / au_2_eV ;
    setm->H_ele_lvcm[3 * 11 + 1] = ( 0.0719 - I * 0.0177 ) / au_2_eV ;
    setm->H_ele_lvcm[4 * 11 + 0] = - ( 0.0719 + I * 0.0177 ) / au_2_eV ;
    setm->H_ele_lvcm[4 * 11 + 2] = ( 0.0719 - I * 0.0177 )  / au_2_eV ;
    setm->H_ele_lvcm[5 * 11 + 1] = - ( 0.0719 + I * 0.0177 ) / au_2_eV ;

    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0086 / au_2_eV    ;
    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.019 / au_2_eV     ;
    setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0086 / au_2_eV    ;
    setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.019 / au_2_eV     ;
    
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0086 / au_2_eV    ;
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.019 / au_2_eV     ;
    setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0086 / au_2_eV    ;
    setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.019 / au_2_eV     ;
     
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0086 / au_2_eV    ;
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.019 / au_2_eV     ;
    setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0086 / au_2_eV     ;
    setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.019 / au_2_eV      ;
    
    //W^{T2, T3}
    setm->H_ele_lvcm[4 * 11 + 8] = ( -0.0270 + I * 0.0046 ) / au_2_eV ;
    setm->H_ele_lvcm[3 * 11 + 9] = - ( -0.0270 - I * 0.0046 ) / au_2_eV ;
    setm->H_ele_lvcm[5 * 11 + 9] = ( -0.0270 + I * 0.0046 ) / au_2_eV ;
    setm->H_ele_lvcm[4 * 11 + 10] = - ( -0.0270 - I * 0.0046 ) / au_2_eV ;
    setm->H_ele_lvcm[8 * 11 + 4] = ( -0.0270 - I * 0.0046 ) / au_2_eV ;
    setm->H_ele_lvcm[9 * 11 + 3] = - ( -0.0270 + I * 0.0046 ) / au_2_eV ;
    setm->H_ele_lvcm[9 * 11 + 5] = ( -0.0270 - I * 0.0046 ) / au_2_eV ;
    setm->H_ele_lvcm[10 * 11 + 4] = - ( -0.0270 + I * 0.0046 ) / au_2_eV ;


}

void sample_bpy(double *P, double *R, struct set_host *setm) {
    int j;
    double x2;

    sample_LVCM(P, R, setm);
    

}

void V_bpy(double *R, double complex *H, int forcetype, struct set_host *setm) {
    
    V_LVCM(R, H, forcetype, setm);
    
}

void dV_bpy(double *R, double complex *dH, int forcetype, struct set_host *setm) {
    
    dV_LVCM(R, dH, forcetype, setm);

}

void nucforce_bpy(double *R, double *nf, struct set_host *setm) {

    nucforce_LVCM(R, nf, setm);

}



