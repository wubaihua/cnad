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


void parameter_pyrazine(double *mass, struct set_host *setm) {
    
    parameter_LVCM(mass, setm);

    switch (setm->type_pyrazine) {
        case 1:

            setm->H_ele_lvcm[0] = 3.94  / au_2_eV ;
            setm->H_ele_lvcm[3] = 4.84  / au_2_eV ;
            setm->omega_lvcm[0] = 0.126 / au_2_eV ;
            setm->omega_lvcm[1] = 0.074 / au_2_eV ;
            setm->omega_lvcm[2] = 0.118 / au_2_eV ;
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.037  / au_2_eV ;
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.105 / au_2_eV ; 
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = -0.254 / au_2_eV ; 
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.149  / au_2_eV ;
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.262  / au_2_eV ;
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.262  / au_2_eV ;
            
            
            break;

        case 2:
            setm->H_ele_lvcm[0] = -0.4617  / au_2_eV ;
            setm->H_ele_lvcm[3] =  0.4617 / au_2_eV ;
            double w[24] = { 0.0936, 0.074,  0.1273, 0.1568, 0.1347,  0.3431, 
                             0.1157, 0.3242, 0.3621, 0.2673, 0.3052,  0.0968,
                             0.0589, 0.04,   0.1726, 0.2863, 0.2484,  0.1536,
                             0.2105, 0.0778, 0.2294, 0.1915, 0.4,     0.381  };
            
            double k1[24] = {   0.0,  -0.0964, 0.047,  0.1594, 0.0308, 0.0782,
                        0.0261, 0.0717, 0.078,  0.056,  0.0625, 0.0188,
                        0.0112, 0.0069, 0.0265, 0.0433, 0.0361, 0.021,
                        0.0281, 0.0102, 0.0284, 0.0196, 0.0306, 0.0269   };

            double k2[24]=  {   0.0,     0.1194, 0.2012,  0.0484,  -0.0308, -0.0782,
                        -0.0261, -0.0717, -0.078,  -0.056,  -0.0625, -0.0188,
                        -0.0112, -0.0069, -0.0265, -0.0433, -0.0361, -0.021,
                        -0.0281, -0.0102, -0.0284, -0.0196, -0.0306, -0.0269   };
            for (int i = 0; i < 24; i++){
                setm->omega_lvcm[i] = w[i] / au_2_eV;
                setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + i] = k1[i] / au_2_eV;
                setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + i] = k2[i] / au_2_eV;
            }
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.1825 / au_2_eV ;
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.1825 / au_2_eV ;
            
     }

     
}

void sample_pyrazine(double *P, double *R, struct set_host *setm) {
    int j;
    double x2;

    sample_LVCM(P, R, setm);
    

}

void V_pyrazine(double *R, double *H, int forcetype, struct set_host *setm) {
    
    V_LVCM(R, H, forcetype, setm);
    
}

void dV_pyrazine(double *R, double *dH, int forcetype, struct set_host *setm) {
    
    dV_LVCM(R, dH, forcetype, setm);

}

void nucforce_pyrazine(double *R, double *nf, struct set_host *setm) {

    nucforce_LVCM(R, nf, setm);

}



