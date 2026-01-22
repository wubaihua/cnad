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


void parameter_soctest(double *mass, struct set_host *setm) {
    
    parameter_LVCM(mass, setm);

    

    
    
    setm->omega_lvcm[0] = 0.01 / au_2_eV ;
    setm->omega_lvcm[1] = 0.02 / au_2_eV ;
    setm->omega_lvcm[2] = 0.03 / au_2_eV ;
    

    //W^{T1, T1}
    setm->H_ele_lvcm[0 * 3 + 0] = 0.0 / au_2_eV ;
    setm->H_ele_lvcm[1 * 3 + 1] = 1.0 / au_2_eV ;
    setm->H_ele_lvcm[2 * 3 + 2] = 2.0 / au_2_eV ;

    setm->H_ele_lvcm[0 * 3 + 1] = (1.0 - I * 0.5) / au_2_eV ;
    setm->H_ele_lvcm[1 * 3 + 0] = (1.0 + I * 0.5) / au_2_eV ;
    setm->H_ele_lvcm[0 * 3 + 2] = (2.0 - I * 1.0) / au_2_eV ;
    setm->H_ele_lvcm[2 * 3 + 0] = (2.0 + I * 1.0) / au_2_eV ;
    setm->H_ele_lvcm[1 * 3 + 2] = (0.1 + I * 2.0) / au_2_eV ;
    setm->H_ele_lvcm[2 * 3 + 1] = (0.1 - I * 2.0) / au_2_eV ;
    

    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.0161  / au_2_eV         ;
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.019  / au_2_eV         ;
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.02 / au_2_eV         ;
    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0002  / au_2_eV         ;
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = -0.0006  / au_2_eV         ;
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.01 / au_2_eV         ;


    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.01  / au_2_eV         ;
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.01  / au_2_eV         ;
    setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.02 / au_2_eV         ;
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.02  / au_2_eV         ;
    setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.03 / au_2_eV         ;
    setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.03 / au_2_eV         ;

    


}

void sample_soctest(double *P, double *R, struct set_host *setm) {
    int j;
    double x2;

    sample_LVCM(P, R, setm);

}

void V_soctest(double *R, double complex *H, int forcetype, struct set_host *setm) {
   
    V_LVCM(R, H, forcetype, setm);
  
    
}

void dV_soctest(double *R, double complex *dH, int forcetype, struct set_host *setm) {
  
    dV_LVCM(R, dH, forcetype, setm);
    
}

void nucforce_soctest(double *R, double *nf, struct set_host *setm) {

    nucforce_LVCM(R, nf, setm);
   

}



void nac_soctest(double *R, double complex *nac, struct set_host *setm) {

    
    memset(nac, 0, setm->Nstate_lvcm * setm->Nstate_lvcm * setm->N_mode_lvcm * sizeof(double complex)); // Initialize nac to zero

    
    for (int k = 0; k < setm->N_mode_lvcm; k++) {
        
        nac[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + k] = 0.2 ;
        nac[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + k] = -0.2 ;
        nac[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + k] = 0.5 ;
        nac[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + k] = -0.5 ;
        nac[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + k] = -0.1 ;
        nac[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + k] = 0.1 ;
        
    }
      
    
}



