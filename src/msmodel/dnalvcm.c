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


void parameter_dnalvcm(double *mass, struct set_host *setm) {

    
    parameter_LVCM(mass, setm);



    switch (setm->type_dnalvcm) {
        case 1: //7HG_CAM

            setm->H_ele_lvcm[0] =    0.00000000e+00   / au_2_eV ;
            setm->H_ele_lvcm[7] =    3.89859110e-01   / au_2_eV ;
            setm->H_ele_lvcm[14] =   6.04562540e-01    / au_2_eV ;
            setm->H_ele_lvcm[21] =   8.20891450e-01    / au_2_eV ;
            setm->H_ele_lvcm[28] =   9.59533130e-01    / au_2_eV ;
            setm->H_ele_lvcm[35] =   1.16304874e+00    / au_2_eV ;


            double wall1[42]  = {4.39558900e-02,  
                            1.79163900e-02,  
                            1.86933900e-02,  
                            2.36805400e-02,  
                            3.83619900e-02,  
                            4.24893800e-02,  
                            4.56578700e-02,  
                            5.10251700e-02,  
                            6.13958200e-02,  
                            6.30969100e-02,  
                            6.78322500e-02,  
                            7.33328300e-02,  
                            7.96819300e-02,  
                            8.07442900e-02,  
                            8.41478400e-02,  
                            9.06885700e-02,  
                            9.27656800e-02,  
                            9.79013100e-02,  
                            1.04880770e-01,  
                            1.10726230e-01,  
                            1.20621690e-01,  
                            1.27405140e-01,  
                            1.36105320e-01,  
                            1.38803210e-01,  
                            1.44309880e-01,  
                            1.49968800e-01,  
                            1.61530860e-01,  
                            1.66971870e-01,  
                            1.73335020e-01,  
                            1.77665000e-01,  
                            1.83557300e-01,  
                            1.87728820e-01,  
                            1.95995750e-01,  
                            1.99568160e-01,  
                            2.02875000e-01,  
                            2.09897640e-01,  
                            2.23841540e-01,  
                            4.04204570e-01,  
                            4.49443230e-01,  
                            4.52616030e-01,  
                            4.56286630e-01,  
                            4.68608650e-01};
        
            for (int i = 0; i < 42; i++){
                setm->omega_lvcm[i] = wall1[i] / au_2_eV ;
            }
        

            //Matrix 0:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.02280703  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = -0.00530551  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = -0.01917402  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.00282707  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.02280703  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.03613622  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.00530551  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.00100436  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.01917402  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = -0.02744961  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.03613622  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.00100436  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = -0.02744961  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.00315864  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.00282707  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.00315864  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            //Matrix 1:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = -0.00280669  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0179683  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.00901538  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.01002689  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.00280669  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.00359302  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0179683  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = -0.00393981  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.00901538  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.00308139  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.00359302  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = -0.00393981  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.00308139  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.00210089  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.01002689  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.00210089  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            //Matrix 2:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = -0.00712585  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = -0.00457975  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = -0.0167397  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.01761203  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.00712585  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = -0.01446392  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.00457975  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.0073124  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.0167397  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = -0.00901327  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = -0.01446392  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0073124  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = -0.00901327  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.01781659  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.01761203  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.01781659  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            //Matrix 3:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.00834819  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = -0.01435411  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = -0.00144752  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = -0.02395962  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.00834819  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0008685  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = -0.01435411  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.00905856  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = -0.00144752  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.01543198  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0008685  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.00905856  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.01543198  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = -0.00549124  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = -0.02395962  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = -0.00549124  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            //Matrix 4:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.01659869  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = -0.01322671  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = -0.01129279  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = -0.00076768  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.00988729  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.00043578  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = -0.00076768  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = -0.0191855  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = -0.00327069  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.01779116  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.00988729  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = -0.00327069  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.02244939  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.00166801  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = -0.01322671  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = -0.01945588  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.00043578  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.01779116  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.00166801  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.01197455  / au_2_eV ;
            
            //Matrix 5:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.02680325  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = -0.00323543  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.02816402  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = -0.00071853  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.00414448  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = -0.00013164  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = -0.00071853  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.00585101  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = -0.00071743  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = -0.0110574  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.00414448  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = -0.00071743  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.01972796  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = -0.00050374  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = -0.00323543  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.02775553  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = -0.00013164  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = -0.0110574  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = -0.00050374  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = -0.01047693  / au_2_eV ;
            
            //Matrix 6:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = -0.01466861  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.00563288  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0218178  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = -0.07709488  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = -0.01466861  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = -0.02021583  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.00563288  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = -0.00781463  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0218178  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.03632237  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = -0.02021583  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = -0.00781463  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.03632237  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.00532711  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = -0.07709488  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.00532711  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            //Matrix 7:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.02174573  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.02199968  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.00931071  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = -0.01469575  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.02174573  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.02743265  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.02199968  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = -0.00768318  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.00931071  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.00100674  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.02743265  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = -0.00768318  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.00100674  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 4.878e-05  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = -0.01469575  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 4.878e-05  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            //Matrix 8:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = -0.01850369  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.00321305  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = -0.03306203  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.00295547  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = -0.00236974  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = -0.00041041  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.00295547  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0497973  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.00529636  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = -0.02586934  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = -0.00236974  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.00529636  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = -0.03510319  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = -0.00290998  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.00321305  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = -0.01510235  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = -0.00041041  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = -0.02586934  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = -0.00290998  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.06816472  / au_2_eV ;
            
            //Matrix 9:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = -0.01255924  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = -0.00192377  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = -0.00104755  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.02674285  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = -0.01255924  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.00628753  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = -0.00192377  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.00144758  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = -0.00104755  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.01123743  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.00628753  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.00144758  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.01123743  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.00062249  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.02674285  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.00062249  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            //Matrix 10:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = -0.05387834  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = -0.00593246  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = -0.05401458  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.00155822  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.00085732  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = -0.00042116  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.00155822  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.00707421  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.00290284  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = -0.01498529  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.00085732  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.00290284  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = -0.04258589  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = -0.00081617  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = -0.00593246  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = -0.0606815  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = -0.00042116  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = -0.01498529  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = -0.00081617  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = -0.02217638  / au_2_eV ;
            
            //Matrix 11:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.05070118  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = -0.01434918  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = -0.00557698  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0144167  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.05070118  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.00130957  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = -0.01434918  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.00132043  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = -0.00557698  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = -0.02294918  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.00130957  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.00132043  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = -0.02294918  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.00080576  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0144167  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.00080576  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            //Matrix 12:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = -0.03347008  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = -0.00843512  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.02775559  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = -0.00082585  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.00303516  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.00030729  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = -0.00082585  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = -0.00435384  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = -0.0020142  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.00147276  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.00303516  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = -0.0020142  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.02163306  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.00044247  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = -0.00843512  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = -0.04830007  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.00030729  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.00147276  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.00044247  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = -0.04557903  / au_2_eV ;
            
            //Matrix 13:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.01737032  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = -0.01260101  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.04848065  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = -0.0187365  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.01737032  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = -0.00675569  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = -0.01260101  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.01017479  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.04848065  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.00246851  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = -0.00675569  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.01017479  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.00246851  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = -0.02398969  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = -0.0187365  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = -0.02398969  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            //Matrix 14:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.01428583  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.01192132  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.03292531  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.00043124  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = -0.01646892  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.00067948  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.00043124  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.01333248  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.00068861  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.02209149  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = -0.01646892  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.00068861  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.00530665  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.00137032  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.01192132  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.00081648  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.00067948  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.02209149  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.00137032  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = -0.03224438  / au_2_eV ;
            
            //Matrix 15:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = -0.03971403  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.02672549  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = -0.03715421  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = -0.00247993  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = -0.03971403  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = -0.07100499  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.02672549  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.01317522  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = -0.03715421  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = -0.03218096  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = -0.07100499  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.01317522  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = -0.03218096  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = -0.02453396  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = -0.00247993  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = -0.02453396  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            //Matrix 16:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = -0.00344527  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.04175581  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = -0.03556431  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = -0.04804924  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = -0.00344527  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.06922991  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.04175581  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = -0.00263618  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = -0.03556431  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.02213407  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.06922991  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = -0.00263618  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.02213407  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = -0.01274088  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = -0.04804924  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = -0.01274088  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            //Matrix 17:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = -0.00861158  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = -0.00648269  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = -0.00708001  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.01832304  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = -0.00861158  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.00264452  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = -0.00648269  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = -0.02523865  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = -0.00708001  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.03363002  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.00264452  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = -0.02523865  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.03363002  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.01083833  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.01832304  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.01083833  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            //Matrix 18:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.00462591  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.00098674  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.00367368  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.00234126  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = -0.00335547  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = -0.00115259  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.00234126  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.07182265  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0060283  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = -0.04397088  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.00335547  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0060283  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = -0.0253072  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = -0.00346201  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.00098674  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = -0.01238115  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.00115259  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = -0.04397088  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = -0.00346201  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.03185342  / au_2_eV ;
            
            //Matrix 19:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.00932217  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = -0.01463371  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.03111497  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.02506949  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.00932217  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = -0.00289372  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = -0.01463371  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.00058972  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.03111497  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = -0.02250206  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = -0.00289372  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.00058972  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = -0.02250206  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.00265795  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.02506949  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.00265795  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            //Matrix 20:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.05782617  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.03246192  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.01986439  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = -0.00146829  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.00696649  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.00134058  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = -0.00146829  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = -0.02365623  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = -0.00565067  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.04890362  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.00696649  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = -0.00565067  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.05591962  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.00419556  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.03246192  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.03469246  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.00134058  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.04890362  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.00419556  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = -0.06178781  / au_2_eV ;
            
            //Matrix 21:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = -0.01918463  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = -0.04144508  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.02108865  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = -0.00044545  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = -0.00327232  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = -0.0007474  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = -0.00044545  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = -0.00748277  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = -0.00102296  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = -0.00862133  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = -0.00327232  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = -0.00102296  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.00911597  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.00052468  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = -0.04144508  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = -0.01047576  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = -0.0007474  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = -0.00862133  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.00052468  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = -0.03306216  / au_2_eV ;
            
            //Matrix 22:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = -0.05075205  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.04367367  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = -0.02748345  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = -0.00080578  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.00046556  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.00104368  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = -0.00080578  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = -0.04638581  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = -0.0024645  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.05254772  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.00046556  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = -0.0024645  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = -0.01333358  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.00503173  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.04367367  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.05047994  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.00104368  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.05254772  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.00503173  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = -0.11225644  / au_2_eV ;
            
            //Matrix 23:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.05210952  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = -0.0280965  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.02312963  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = -0.00044414  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.00183922  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = -0.00012382  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = -0.00044414  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.01945682  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.00119519  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.00967884  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.00183922  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.00119519  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.00108843  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = -0.00148872  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = -0.0280965  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.06108975  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = -0.00012382  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.00967884  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = -0.00148872  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.00149598  / au_2_eV ;
            
            //Matrix 24:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = -0.0425859  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = -0.01116436  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.02666737  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = -0.00099204  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = -0.00758355  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.00015865  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = -0.00099204  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = -0.03047434  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = -0.00120384  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.02767796  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = -0.00758355  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = -0.00120384  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = -0.01197319  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0030246  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = -0.01116436  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = -0.00217679  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.00015865  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.02767796  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0030246  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = -0.08653456  / au_2_eV ;
            
            //Matrix 25:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = -0.09007095  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.03723079  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = -0.03564648  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.00179093  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = -0.0136204  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = -0.00090169  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.00179093  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.00122317  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.004764  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = -0.03154589  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = -0.0136204  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.004764  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = -0.06163452  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = -0.00429516  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.03723079  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = -0.04054359  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = -0.00090169  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = -0.03154589  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = -0.00429516  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.01986601  / au_2_eV ;
            
            //Matrix 26:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.07796021  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.02136097  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0455792  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = -0.00151506  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = -0.00263905  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = -0.0001088  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = -0.00151506  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = -0.02897984  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = -0.0027347  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.00525018  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = -0.00263905  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = -0.0027347  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.02095282  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.00054892  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.02136097  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.01170123  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = -0.0001088  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.00525018  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.00054892  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.00476149  / au_2_eV ;
            
            //Matrix 27:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = -0.01469404  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.01077076  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.00231239  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = -0.00111054  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = -0.00522251  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0001889  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = -0.00111054  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = -0.0406817  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = -0.00212879  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = -0.00940949  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = -0.00522251  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = -0.00212879  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = -0.00408101  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = -0.00070966  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.01077076  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.00217681  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0001889  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = -0.00940949  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = -0.00070966  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = -0.034966  / au_2_eV ;
            
            //Matrix 28:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.00612194  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0349676  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.017279  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.00071353  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = -0.00724514  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = -0.00059352  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.00071353  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = -0.01973248  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.00302711  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = -0.03113125  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.00724514  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.00302711  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = -0.05741629  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = -0.00277939  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0349676  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.05605602  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.00059352  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = -0.03113125  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = -0.00277939  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.02707965  / au_2_eV ;
            
            //Matrix 29:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = -0.17564413  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.07238678  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = -0.03197297  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = -0.00128566  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.01300904  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.00247524  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = -0.00128566  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = -0.09181221  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.0007417  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0810006  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = -0.01300904  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = -0.0007417  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.08898128  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.00506208  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.07238678  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = -0.09293201  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.00247524  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0810006  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.00506208  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = -0.16356667  / au_2_eV ;
            
            //Matrix 30:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = -0.0035376  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = -0.0156158  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.00734751  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.00033252  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = -0.00779635  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.00121583  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.00033252  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.01959201  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.00085367  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.04919609  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = -0.00779635  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.00085367  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0197276  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.00500442  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = -0.0156158  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = -0.05578317  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.00121583  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.04919609  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.00500442  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.01415028  / au_2_eV ;
            
            //Matrix 31:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.09904899  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = -0.00852294  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.08952547  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = -0.00113738  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.01593388  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = -0.00019155  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = -0.00113738  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.05714373  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = -0.00165517  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.00081506  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.01593388  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = -0.00165517  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.07836868  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.00210812  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = -0.00852294  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.1025872  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = -0.00019155  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.00081506  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.00210812  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.01306154  / au_2_eV ;
            
            //Matrix 32:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = -0.00367324  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.01158926  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = -0.00979563  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.00012874  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.02977802  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.00017332  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.00012874  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.05864134  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.00035935  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.00654381  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.02977802  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.00035935  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.04095265  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.00240415  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.01158926  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.08204196  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.00017332  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.00654381  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.00240415  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = -0.01483105  / au_2_eV ;
            
            //Matrix 33:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.12435135  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.09718555  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = 0.13170659  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = -0.00398248  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.01346995  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = -0.0015611  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = -0.00398248  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.01996939  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = -0.00762887  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = -0.09567214  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = 0.01346995  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = -0.00762887  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.14598594  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = -0.00608694  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.09718555  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.16068766  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = -0.0015611  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = -0.09567214  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = -0.00608694  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = 0.07608613  / au_2_eV ;
            
            //Matrix 34:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.01537503  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.02475694  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = -0.06680268  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.00040921  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0052102  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = -0.00139164  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.00040921  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = -0.01387831  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.00155001  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = -0.02149587  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.0052102  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.00155001  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = -0.03646448  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = -0.00269829  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.02475694  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = -0.11360804  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = -0.00139164  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = -0.02149587  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = -0.00269829  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = 0.05809682  / au_2_eV ;
            
            //Matrix 35:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.06816512  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = -0.01296792  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = 0.19319803  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = -0.0025803  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = -0.02348189  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.00010407  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = -0.0025803  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.03972787  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = -0.00414496  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.00068979  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = -0.02348189  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = -0.00414496  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.12191018  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.0012778  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = -0.01296792  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.2039485  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = 0.00010407  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.00068979  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0012778  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.01374183  / au_2_eV ;
            
            //Matrix 36:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.08598742  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.00996145  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = -0.10083527  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.01338015  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = -0.00502317  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.00415137  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.01338015  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.32136857  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.02444536  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.13598764  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = -0.00502317  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.02444536  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = -0.10155231  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.00395807  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.00996145  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = -0.04829968  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.00415137  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.13598764  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.00395807  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.21177312  / au_2_eV ;
            
            //Matrix 37:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0149662  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = -0.01814611  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.00060916  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0013867  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = -0.00057687  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.00060916  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.02367336  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.00139293  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = -0.01474413  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.0013867  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.00139293  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.007075  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = -0.00252885  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = -0.01814611  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.01605475  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = -0.00057687  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = -0.01474413  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = -0.00252885  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0304772  / au_2_eV ;
            
            //Matrix 38:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = -0.01673507  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = -0.00968993  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.05726812  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.0010199  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = -0.02489567  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.00275004  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.0010199  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = -0.02340176  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = -0.00191889  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = -0.00223752  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = -0.02489567  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = -0.00191889  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.01089621  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = -0.00105033  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = -0.00968993  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = -0.0035374  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.00275004  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = -0.00223752  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = -0.00105033  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = -0.01088442  / au_2_eV ;
            
            //Matrix 39:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 39] = 0.01564652  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 39] = -0.00636692  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 39] = 0.12393472  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 39] = -0.00115344  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 39] = -0.0503214  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 39] = -0.00013876  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 39] = -0.00115344  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 39] = -0.00340149  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 39] = -0.00116665  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 39] = -0.00187463  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 39] = -0.0503214  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 39] = -0.00116665  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 39] = 0.04559213  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 39] = 0.00026015  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 39] = -0.00636692  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 39] = 0.00639469  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 39] = -0.00013876  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 39] = -0.00187463  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 39] = 0.00026015  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 39] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 39] = -0.00408169  / au_2_eV ;
            
            //Matrix 40:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 40] = 0.00897975  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 40] = 0.00111026  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 40] = -0.02203486  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 40] = 0.0021022  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 40] = -0.01454676  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 40] = -0.00089341  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 40] = 0.0021022  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 40] = 0.01415622  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 40] = 0.0101406  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 40] = -0.00325214  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 40] = -0.01454676  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 40] = 0.0101406  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 40] = -0.13593282  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 40] = -0.00686132  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 40] = 0.00111026  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 40] = -0.01006821  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 40] = -0.00089341  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 40] = -0.00325214  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 40] = -0.00686132  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 40] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 40] = 0.01945553  / au_2_eV ;
            
            //Matrix 41:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 41] = -0.00204086  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 41] = 0.00096732  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 41] = -0.0508702  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 41] = 0.00016319  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 41] = 0.0284932  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 41] = -6.815e-05  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 41] = 0.00016319  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 41] = 0.00149663  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 41] = 0.00022216  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 41] = 0.000956  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 41] = 0.0284932  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 41] = 0.00022216  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 41] = -0.01960723  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 41] = 7.002e-05  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 41] = 0.00096732  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 41] = -0.00176874  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 41] = -6.815e-05  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 41] = 0.000956  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 41] = 7.002e-05  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 41] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 41] = -0.0009524  / au_2_eV ;





            break;

        case 2: //C_CAM
            setm->H_ele_lvcm[0] =     0.00000000e+00   / au_2_eV ;
            setm->H_ele_lvcm[8] =     2.72579460e-01   / au_2_eV ;
            setm->H_ele_lvcm[16] =    7.95362710e-01    / au_2_eV ;
            setm->H_ele_lvcm[24] =    8.98675120e-01    / au_2_eV ;
            setm->H_ele_lvcm[32] =    9.25141730e-01    / au_2_eV ;
            setm->H_ele_lvcm[40] =    1.11982099e+00    / au_2_eV ;
            setm->H_ele_lvcm[48] =    1.32986639e+00    / au_2_eV ;                            

            double wall2[33]  = {   
                                9.39475000e-03,  
                                1.64171900e-02,  
                                2.49513800e-02,  
                                4.50746700e-02,  
                                5.02487700e-02,  
                                6.69148500e-02,  
                                6.70268100e-02,  
                                6.81703400e-02,  
                                7.25531500e-02,  
                                7.81679600e-02,  
                                9.10407400e-02,  
                                9.70918400e-02,  
                                9.72133100e-02,  
                                9.87367200e-02,  
                                1.16450630e-01,  
                                1.22601560e-01,  
                                1.23822680e-01,  
                                1.36787540e-01,  
                                1.41166690e-01,  
                                1.51984490e-01,  
                                1.59266320e-01,  
                                1.70025960e-01,  
                                1.81227580e-01,  
                                1.88795180e-01,  
                                1.98542320e-01,  
                                2.03679860e-01,  
                                2.13507290e-01,  
                                2.24021320e-01,  
                                3.98869130e-01,  
                                4.02112550e-01,  
                                4.50410650e-01,  
                                4.52498940e-01,  
                                4.67204480e-01};                      
                                
        
            for (int i = 0; i < 33; i++){
                setm->omega_lvcm[i] = wall2[i] / au_2_eV ;
            }


            //Matrix 0:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = -0.00046166  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = -0.01087625  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = -0.0072752  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = -0.00374822  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = -0.0239755  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.00046166  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.10540552  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.01087625  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.01540588  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.0072752  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = -0.02485543  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.10540552  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.01540588  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = -0.02485543  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.03102959  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = -0.04925849  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.00374822  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.03102959  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.0239755  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = -0.04925849  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;
            
            //Matrix 1:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = -0.00132338  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.01374157  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.04400769  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = -0.04802835  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.01122644  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.00132338  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.0125695  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.01374157  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = -0.00534115  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.04400769  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = -0.04957658  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0125695  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = -0.00534115  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = -0.04957658  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.03476724  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.00211596  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.04802835  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.03476724  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.01122644  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.00211596  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;
            
            //Matrix 2:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = -0.05251693  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = -0.01583769  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = -0.04149005  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.03285914  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = -0.00706581  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.05251693  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = -0.05463118  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.01583769  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.00963998  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.04149005  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = -0.0281166  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = -0.05463118  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.00963998  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = -0.0281166  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.02505558  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.00428379  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.03285914  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.02505558  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.00706581  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.00428379  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;
            
            //Matrix 3:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = -0.013741  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = -0.03650207  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = -0.03687085  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.00352332  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = -0.02706896  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.01148263  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = -0.00089624  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.00352332  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.01850604  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.00766441  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = -0.00494031  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = -0.02458275  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = -0.02706896  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.00766441  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = -0.03825334  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.04180983  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = -0.00495534  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = -0.03650207  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = -0.05496769  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.01148263  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = -0.00494031  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.04180983  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = -0.00773897  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.0009339  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = -0.00089624  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = -0.02458275  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = -0.00495534  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0009339  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = -0.00054197  / au_2_eV ;
            
            //Matrix 4:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0319297  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = -0.05295958  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = -0.05967554  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.04432138  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.00354583  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0319297  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = -0.0183656  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = -0.05295958  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.00205785  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = -0.05967554  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.01106529  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = -0.0183656  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.00205785  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.01106529  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = -5.257e-05  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = -0.04747372  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.04432138  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = -5.257e-05  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.00354583  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = -0.04747372  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;
            
            //Matrix 5:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = -0.04322141  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.00600601  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = -0.0031452  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.00216434  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.00149728  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = -0.04322141  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.02200817  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.00600601  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = -0.00106137  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = -0.0031452  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.02442674  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.02200817  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = -0.00106137  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.02442674  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = -0.01672436  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = -0.01903857  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.00216434  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = -0.01672436  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.00149728  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = -0.01903857  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;
            
            //Matrix 6:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.06993317  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = -0.00350528  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.02789067  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.00415557  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = -0.03479185  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.02557977  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = -0.00062931  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.00415557  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.02857263  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.00245143  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.00032105  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.01462553  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = -0.03479185  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.00245143  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.00734854  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = -0.0012311  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.00143734  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = -0.00350528  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.01156486  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.02557977  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.00032105  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = -0.0012311  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.04598609  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = -0.0007457  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = -0.00062931  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.01462553  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.00143734  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = -0.0007457  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = -0.0131976  / au_2_eV ;
            
            //Matrix 7:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = -0.01455858  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = -0.02230618  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = -0.07333146  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0057143  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = -0.04907968  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.04442815  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = -0.00454992  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0057143  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = -0.03020387  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0064763  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = -0.00160327  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = -0.02315703  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = -0.04907968  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0064763  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = -0.07551769  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.00845257  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = -0.00581567  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = -0.02230618  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.05619197  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.04442815  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = -0.00160327  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.00845257  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = -0.00693772  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.00140797  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = -0.00454992  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = -0.02315703  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = -0.00581567  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.00140797  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.0254437  / au_2_eV ;
            
            //Matrix 8:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = -0.04952506  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.01195284  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.15661022  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = -0.00577572  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.03366883  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = -0.03406625  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.00396071  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = -0.00577572  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = -0.04707594  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = -0.0015196  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = -0.00319676  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.00401632  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.03366883  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = -0.0015196  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = -0.02913429  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.03214236  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.00023007  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.01195284  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.00639505  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = -0.03406625  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = -0.00319676  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.03214236  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = -0.02720162  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.00027663  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.00396071  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.00401632  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.00023007  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.00027663  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = -0.00217706  / au_2_eV ;
            
            //Matrix 9:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0731667  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.03989163  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.03376095  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = -0.00747654  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.01054599  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0731667  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = -0.05530502  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.03989163  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = -0.01353924  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.03376095  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = -0.00249299  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = -0.05530502  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = -0.01353924  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = -0.00249299  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = -0.00933853  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.02547945  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = -0.00747654  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = -0.00933853  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.01054599  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.02547945  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;
            
            //Matrix 10:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = -0.04948964  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.00739549  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = -0.02364165  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.00941533  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0024672  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = -0.04948964  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.03252363  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.00739549  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.01422633  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = -0.02364165  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = -0.0186987  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.03252363  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.01422633  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = -0.0186987  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.02048872  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = -0.00829687  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.00941533  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.02048872  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0024672  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = -0.00829687  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;
            
            //Matrix 11:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0457152  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = -0.04114897  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.06231382  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 6.322e-05  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = -0.00333186  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = -0.0192147  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.00134518  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 6.322e-05  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = -0.01142801  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = -1.378e-05  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = -0.00525517  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = -0.0165716  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = -0.00333186  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = -1.378e-05  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = -0.00025388  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.04559575  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = -0.00261738  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = -0.04114897  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.06462687  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = -0.0192147  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = -0.00525517  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.04559575  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = -0.02913422  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.00137542  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.00134518  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = -0.0165716  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = -0.00261738  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.00137542  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.00476125  / au_2_eV ;
            
            //Matrix 12:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.01613423  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = -0.01455818  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.00201106  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.01297696  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = -0.0105487  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.01613423  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.014387  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = -0.01455818  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.01153127  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.00201106  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = -0.0174621  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.014387  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.01153127  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = -0.0174621  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.01485686  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = -0.01749842  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.01297696  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.01485686  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = -0.0105487  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = -0.01749842  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;
            
            //Matrix 13:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.03357166  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = -0.00450921  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.02137594  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.00635584  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.00164564  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.03357166  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.04350052  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = -0.00450921  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0001027  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.02137594  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.02234387  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.04350052  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0001027  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.02234387  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = -0.00957026  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.00462293  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.00635584  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = -0.00957026  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.00164564  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.00462293  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;
            
            //Matrix 14:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.03415007  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.00339913  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = -0.04598862  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0064121  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = -0.04660165  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.02563984  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = -0.00264648  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0064121  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.00870897  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.00538113  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = -0.00254098  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = -0.00559284  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = -0.04660165  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.00538113  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = -0.03102367  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.0268257  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = -0.00289548  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.00339913  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.02979663  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.02563984  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = -0.00254098  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0268257  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = -0.0272086  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.00174354  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = -0.00264648  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = -0.00559284  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = -0.00289548  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.00174354  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.00952402  / au_2_eV ;
            
            //Matrix 15:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.04579311  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = -0.00165552  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = -7.561e-05  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.00597751  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 4.036e-05  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.04579311  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = -0.01760854  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = -0.00165552  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.00570124  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = -7.561e-05  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = -0.00104303  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = -0.01760854  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.00570124  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = -0.00104303  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0040626  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.01529514  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.00597751  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0040626  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 4.036e-05  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.01529514  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;
            
            //Matrix 16:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = -0.04666779  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.01535227  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = -0.13905687  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.00617623  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = -0.04592395  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.02040554  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = -0.00430231  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.00617623  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = -0.00394511  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.00194256  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = -0.00434619  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.00998619  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = -0.04592395  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.00194256  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = -0.01756456  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.04651091  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = -0.00112024  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.01535227  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = -0.01836736  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.02040554  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = -0.00434619  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.04651091  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.02219699  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.00214069  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = -0.00430231  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.00998619  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = -0.00112024  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.00214069  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.00965981  / au_2_eV ;
            
            //Matrix 17:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = -0.02122528  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.04539314  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = -0.06000035  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = -0.00114796  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.01137038  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.02771762  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = -0.0032701  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = -0.00114796  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = -0.02339588  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = -0.00197047  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.00266967  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.03276302  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.01137038  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = -0.00197047  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = -0.01865218  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = -0.02001584  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.00314273  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.04539314  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0435386  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.02771762  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.00266967  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = -0.02001584  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.06899268  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.00023853  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = -0.0032701  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.03276302  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.00314273  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.00023853  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.00938185  / au_2_eV ;
            
            //Matrix 18:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = -0.08081818  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = -0.01375444  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.09306402  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0008429  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = -0.00309018  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = -0.0325584  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = -0.00117879  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0008429  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.04410097  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.01280462  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = -0.00744182  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.01089797  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.00309018  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.01280462  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = -0.05468888  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.05272229  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = -0.0007762  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = -0.01375444  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = -0.0874841  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.0325584  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = -0.00744182  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.05272229  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = -0.08805201  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.00211368  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.00117879  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.01089797  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = -0.0007762  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.00211368  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = -0.0078914  / au_2_eV ;
            
            //Matrix 19:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.02979593  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = -0.02794765  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = -0.02639505  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.00220926  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = -0.00262592  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = -0.00261107  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0009344  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.00220926  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.09076661  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.01380487  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = -0.00191179  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.00331655  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = -0.00262592  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.01380487  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = -0.02505354  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.01532652  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = -0.00056878  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = -0.02794765  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = -0.00870713  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = -0.00261107  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = -0.00191179  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.01532652  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = -0.0048967  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 9.977e-05  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0009344  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.00331655  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = -0.00056878  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 9.977e-05  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.03034171  / au_2_eV ;
            
            //Matrix 20:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.1239508  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = -0.065199  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.05361729  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = -7.979e-05  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = -0.00728331  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = -0.0462687  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.00452774  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = -7.979e-05  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0254256  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = -0.01016786  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.01023446  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = -0.03549332  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = -0.00728331  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = -0.01016786  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.11447036  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = -0.10309197  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.00163352  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = -0.065199  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.00816036  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = -0.0462687  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.01023446  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = -0.10309197  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = -0.10004305  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = -0.00336774  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.00452774  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = -0.03549332  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.00163352  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = -0.00336774  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = -0.01360484  / au_2_eV ;
            
            //Matrix 21:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.07115828  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = -0.02009656  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0151028  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = -0.001906  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.01852065  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = -0.00396963  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.00181536  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = -0.001906  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.06244832  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.00394647  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.00119019  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = -0.0111301  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.01852065  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.00394647  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.02775462  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = -0.01244395  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.0002291  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = -0.02009656  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = -0.00911638  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = -0.00396963  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.00119019  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = -0.01244395  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.04476304  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = -0.00182528  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.00181536  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = -0.0111301  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0002291  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = -0.00182528  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = -0.01169901  / au_2_eV ;
            
            //Matrix 22:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = -0.08149826  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = -0.04913099  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = -0.03319743  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.00144177  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = -0.00178983  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = -0.03400067  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.001273  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.00144177  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.01047882  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.00367857  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = -0.00302907  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = -0.01730225  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = -0.00178983  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.00367857  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = -0.01346267  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.020372  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = -0.00165419  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = -0.04913099  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = -0.12639654  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = -0.03400067  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = -0.00302907  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.020372  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = -0.10300205  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.00088899  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.001273  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = -0.01730225  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = -0.00165419  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.00088899  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = -0.03660207  / au_2_eV ;
            
            //Matrix 23:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = -0.04408258  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = -0.02171413  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = -0.02898192  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.00756113  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = -0.05966256  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.02740394  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = -0.00236583  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.00756113  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = -0.0749684  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = -0.00661227  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = -0.00207084  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = -0.00497316  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = -0.05966256  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = -0.00661227  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = -0.01591928  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.01626709  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = -0.00200633  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = -0.02171413  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.00952415  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.02740394  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = -0.00207084  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.01626709  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = -0.05809357  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.00197878  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = -0.00236583  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = -0.00497316  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = -0.00200633  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.00197878  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = -0.0312922  / au_2_eV ;
            
            //Matrix 24:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.26693976  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.08772773  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.17187284  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = -0.00150625  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.01457163  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.08919161  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.00070497  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = -0.00150625  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.01789533  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = -0.01027995  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.02619653  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.02937392  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.01457163  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = -0.01027995  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.08444566  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = -0.21971747  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.00815172  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.08772773  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.29592726  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.08919161  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.02619653  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = -0.21971747  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.17993679  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = -0.00548779  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.00070497  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.02937392  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.00815172  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = -0.00548779  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.04300096  / au_2_eV ;
            
            //Matrix 25:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = -0.05102095  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = -0.0108431  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = -0.01156537  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.00255426  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = -0.02536383  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = -0.00227832  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = -0.00052951  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.00255426  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = -0.05347145  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = -0.00303671  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = -0.00129898  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = -0.01028033  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = -0.02536383  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = -0.00303671  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = -0.02680204  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.01062489  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = -0.00244059  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = -0.0108431  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = -0.0048984  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = -0.00227832  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = -0.00129898  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.01062489  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = -0.07755329  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.00115123  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = -0.00052951  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = -0.01028033  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = -0.00244059  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.00115123  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = -0.02054316  / au_2_eV ;
            
            //Matrix 26:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.00680182  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = -0.04970917  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = -0.14000755  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.00128666  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.00485426  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.04007585  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = -0.00251154  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.00128666  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.07513802  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.02711723  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = -0.01235916  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = -0.02868969  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.00485426  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.02711723  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = -0.14141411  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.09795599  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = -0.00906706  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = -0.04970917  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.02707632  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.04007585  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = -0.01235916  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.09795599  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.06506005  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.00229246  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = -0.00251154  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = -0.02868969  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = -0.00906706  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.00229246  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.0236705  / au_2_eV ;
            
            //Matrix 27:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = -0.15264919  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = -0.04019132  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.06403411  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = -0.00864692  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0740705  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = -0.00567574  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.00113038  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = -0.00864692  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = -0.01498328  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0563041  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = -0.00444455  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = -0.02964177  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0740705  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0563041  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = -0.50036289  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.01061566  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = -0.01622682  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = -0.04019132  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = -0.02980299  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = -0.00567574  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = -0.00444455  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.01061566  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = -0.36815772  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.00675573  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.00113038  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = -0.02964177  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = -0.01622682  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.00675573  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.04000031  / au_2_eV ;
            
            //Matrix 28:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = -0.03850407  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = -0.00033436  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.02598686  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.00079852  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = -0.0027881  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = -0.00938173  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = -0.00048665  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.00079852  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.02149601  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.00558667  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = -0.0030807  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.00874236  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.0027881  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.00558667  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = -0.01482725  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.02165498  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0003121  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = -0.00033436  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = -0.03904831  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.00938173  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = -0.0030807  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.02165498  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = -0.02422062  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0009824  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.00048665  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.00874236  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0003121  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0009824  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = -0.00326482  / au_2_eV ;
            
            //Matrix 29:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.01102063  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = -0.00190827  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = -0.00231293  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.00032617  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.00204799  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.00521309  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = -0.00085767  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.00032617  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = -0.04040737  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.00388065  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = -0.00042039  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.00491264  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.00204799  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = -0.00388065  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.00544295  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.00251666  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = -0.00018953  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = -0.00190827  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0156465  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.00521309  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = -0.00042039  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.00251666  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.02231331  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.00048842  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = -0.00085767  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.00491264  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.00018953  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.00048842  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.00680198  / au_2_eV ;
            
            //Matrix 30:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.01537446  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = -0.00419633  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.02544284  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = -0.00018305  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = -0.00016607  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = -0.00448409  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.00672162  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = -0.00018305  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = -0.0097948  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = -0.0044327  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.00082427  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.02236336  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = -0.00016607  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = -0.0044327  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.02204228  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = -0.00615764  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.00744793  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = -0.00419633  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = -0.03823201  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = -0.00448409  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.00082427  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = -0.00615764  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.00258591  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = -0.00377993  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.00672162  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.02236336  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.00744793  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = -0.00377993  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = -0.15864566  / au_2_eV ;
            
            //Matrix 31:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = -0.00013584  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = -0.01484373  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = -0.02312961  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.00163411  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.01167523  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = -0.00027663  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = -6.945e-05  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.00163411  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.12831217  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0144591  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = -0.00115227  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.02525164  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.01167523  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0144591  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = -0.03091105  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.00396214  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 3.046e-05  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = -0.01484373  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = -0.00217712  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = -0.00027663  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = -0.00115227  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.00396214  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = -0.00421777  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.00046219  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = -6.945e-05  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.02525164  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 3.046e-05  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.00046219  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.04709119  / au_2_eV ;
            
            //Matrix 32:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = -0.00054421  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = -0.00348874  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = -0.00312928  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 1.912e-05  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = -0.00118704  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = -0.00118729  / au_2_eV ;
            
            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.00373132  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 1.912e-05  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = -0.01441993  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = -0.00187305  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = -0.00043005  / au_2_eV ;
            
            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.01385295  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = -0.00118704  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = -0.00187305  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 2.8e-07  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.00146524  / au_2_eV ;
            
            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.00477742  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = -0.00348874  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = -0.00625863  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = -0.00118729  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = -0.00043005  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.00146524  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = -0.00353732  / au_2_eV ;
            
            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = -0.00171015  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.00373132  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.01385295  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.00477742  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = -0.00171015  / au_2_eV ;
            
            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = -0.08136451  / au_2_eV ;



            break;
        case 3: //T_CAM

            setm->H_ele_lvcm[0]   =   0.00000000e+00     / au_2_eV ;
            setm->H_ele_lvcm[8]   =   1.70946170e-01     / au_2_eV ;
            setm->H_ele_lvcm[16]  =   8.02312560e-01      / au_2_eV ;
            setm->H_ele_lvcm[24]  =   1.32490578e+00      / au_2_eV ;
            setm->H_ele_lvcm[32]  =   1.53364164e+00      / au_2_eV ;
            setm->H_ele_lvcm[40]  =   1.58867100e+00      / au_2_eV ;
            setm->H_ele_lvcm[48]  =   1.63653790e+00      / au_2_eV ;

            double wall3[39]  = {
                                1.34386100e-02,  
                                1.81761000e-02,  
                                1.86040900e-02,  
                                3.47238900e-02,  
                                3.66405500e-02,  
                                4.92234700e-02,  
                                4.95552800e-02,  
                                5.78266100e-02,  
                                6.87386800e-02,  
                                6.88042200e-02,  
                                7.63024000e-02,  
                                8.31562200e-02,  
                                9.23884500e-02,  
                                9.59022000e-02,  
                                9.67055900e-02,  
                                1.01228630e-01,  
                                1.15128920e-01,  
                                1.21583470e-01,  
                                1.28157130e-01,  
                                1.33774230e-01,  
                                1.45426940e-01,  
                                1.50864030e-01,  
                                1.54952340e-01,  
                                1.72426720e-01,  
                                1.76595940e-01,  
                                1.77950380e-01,  
                                1.79616710e-01,  
                                1.83088290e-01,  
                                1.86172200e-01,  
                                1.89064230e-01,  
                                2.15767460e-01,  
                                2.22741850e-01,  
                                2.27682310e-01,  
                                3.79518740e-01,  
                                3.86914520e-01,  
                                3.89200090e-01,  
                                3.98980800e-01,  
                                4.49365730e-01,  
                                4.54687430e-01};
        
            for (int i = 0; i < 39; i++){
                setm->omega_lvcm[i] = wall3[i] / au_2_eV ;
            }

            //Matrix 0:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.00958436  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = -0.00984629  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = -0.00600596  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.00958436  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = -0.01005279  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.01943114  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = -0.01970208  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = -0.01005279  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.03713601  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.01943114  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = -0.01141918  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.02290332  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.00984629  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = -0.01141918  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.00152331  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = -0.00600596  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.03713601  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.02290332  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = -0.0042922  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 0] = -0.01970208  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 0] = 0.00152331  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 0] = -0.0042922  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 0] = 0.0  / au_2_eV ;

            //Matrix 1:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.01369512  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = -0.02831022  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.00385648  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.01369512  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = -0.00869794  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.00307896  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.00254946  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = -0.00869794  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = -0.00263599  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = -0.00446663  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.00307896  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.01349455  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = -0.00472465  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = -0.02831022  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = -0.00263599  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.01349455  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.00675049  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.00385648  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = -0.00446663  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = -0.00472465  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.02727938  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 1] = 0.00254946  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 1] = 0.00675049  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 1] = 0.02727938  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 1] = 0.0  / au_2_eV ;

            //Matrix 2:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = -0.01755137  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = -0.0346372  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.00946139  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.01755137  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = -0.02395479  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.00893727  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = -0.02671857  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = -0.02395479  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.00145582  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.00492434  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.00893727  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.00693447  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = -0.00985724  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = -0.0346372  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.00145582  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.00693447  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.001062  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.00946139  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.00492434  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = -0.00985724  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = -0.01424537  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 2] = -0.02671857  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 2] = 0.001062  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 2] = -0.01424537  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 2] = 0.0  / au_2_eV ;

            //Matrix 3:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = -0.00666679  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.00380958  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = -0.00235254  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.00598625  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.00596753  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.00517016  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.01224507  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = -0.00238202  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = -0.00235254  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = -0.00238202  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = -0.00149657  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 3] = 0.00596753  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 3] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 3] = -0.01265304  / au_2_eV ;

            //Matrix 4:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = -0.00944081  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = -0.03192161  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = -0.00944081  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = -0.00980859  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.02013987  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = -0.0450863  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = -0.00980859  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.00345829  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.02719198  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.02013987  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = -0.00910539  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.02191833  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = -0.03192161  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.00345829  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = -0.00910539  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.02719198  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.02191833  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = -0.04394144  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 4] = -0.0450863  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 4] = -0.04394144  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 4] = 0.0  / au_2_eV ;

            //Matrix 5:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = -0.01510235  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = -0.00365297  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.03061275  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = -0.00522461  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.02911611  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = -0.00498711  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = -0.00365297  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.00244905  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = -0.00522461  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = -0.06520498  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.00925521  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.00925521  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.05010271  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 5] = -0.00498711  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 5] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 5] = 0.03918442  / au_2_eV ;

            //Matrix 6:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = -0.00331285  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.02800787  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = -0.00331285  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = -0.06722531  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.02060804  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.11506202  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = -0.06722531  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = -0.00177189  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = -0.059144  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.02060804  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.00767123  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.04855619  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.02800787  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = -0.00177189  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.00767123  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = -0.00694155  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = -0.059144  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.04855619  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = -0.0360073  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 6] = 0.11506202  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 6] = -0.00694155  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 6] = -0.0360073  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 6] = 0.0  / au_2_eV ;

            //Matrix 7:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.00448987  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = -0.03102213  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = -0.04080917  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.00980009  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = -0.03782384  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = -0.00781451  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = -0.01537442  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = -0.04080917  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.00612653  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = -0.00343855  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.00980009  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = -0.00343855  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = -0.03047952  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 7] = -0.00781451  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 7] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 7] = -0.02462624  / au_2_eV ;

            //Matrix 8:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = -0.00598644  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = -0.00495699  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = -0.04489891  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.01773027  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = -0.00749553  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = -0.03061272  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.01182625  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = -0.00495699  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = -0.00612261  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.01773027  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.01945677  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.00121618  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = -0.00749553  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.00121618  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = -0.0172797  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 8] = 0.01182625  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 8] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 8] = -0.05292616  / au_2_eV ;

            //Matrix 9:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = -0.0067466  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = -0.01769227  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = -0.0067466  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.00571838  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.00238413  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = -0.01782152  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.00571838  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = -0.07151789  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.00238413  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = -0.00474137  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = -0.00548745  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = -0.01769227  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = -0.00474137  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0016445  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = -0.07151789  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = -0.00548745  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = -0.03198536  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 9] = -0.01782152  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 9] = 0.0016445  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 9] = -0.03198536  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 9] = 0.0  / au_2_eV ;

            //Matrix 10:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0160547  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = -0.03986491  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = -0.02048946  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.02615677  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = -0.00775525  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0079557  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = -0.02048946  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = -0.0435383  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.00126514  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.02615677  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.00126514  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = -0.01646249  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 10] = 0.0079557  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 10] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 10] = -0.01877583  / au_2_eV ;

            //Matrix 11:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.01558802  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = -0.00836029  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0030363  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = -0.00975059  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = -0.00374039  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0030363  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = -0.00233614  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = -0.01330198  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = -0.00975059  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = -0.00310946  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.00163476  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.01558802  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = -0.00233614  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = -0.00310946  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = -0.00836029  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = -0.01330198  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.00163476  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.00629059  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 11] = -0.00374039  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 11] = 0.00629059  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 11] = 0.0  / au_2_eV ;

            //Matrix 12:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = -0.01578267  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.00862478  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = -0.04639572  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = -0.00918769  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.00394538  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = -0.01278934  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.00862478  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = -0.04870826  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = -0.00918769  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = -0.01462149  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.00394538  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = -0.01462149  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = -0.08645605  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 12] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 12] = -0.00734707  / au_2_eV ;

            //Matrix 13:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.01485602  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0420459  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = -0.01039101  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.01485602  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.00934065  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.008768  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = -0.00638581  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.00934065  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.00643571  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0764214  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.008768  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.00873602  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = -0.01836902  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0420459  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.00643571  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.00873602  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = -0.01238362  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = -0.01039101  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0764214  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = -0.01836902  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = -0.13272137  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 13] = -0.00638581  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 13] = -0.01238362  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 13] = -0.13272137  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 13] = 0.0  / au_2_eV ;

            //Matrix 14:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = -0.01317936  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = -0.01501286  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.02981283  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = -0.01317936  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = -0.0136652  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0083795  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = -0.0440355  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = -0.0136652  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = -0.02548263  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0083795  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = -0.00556295  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = -0.00128774  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = -0.01501286  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = -0.00556295  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.01111411  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.02981283  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = -0.02548263  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = -0.00128774  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.09043536  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 14] = -0.0440355  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 14] = 0.01111411  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 14] = 0.09043536  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 14] = 0.0  / au_2_eV ;

            //Matrix 15:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = -0.00340147  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = -0.00697675  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.05932043  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = -0.01749499  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = -0.01176863  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.00435382  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.00241213  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = -0.00697675  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.00857163  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = -0.01749499  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.04404997  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.01659885  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = -0.01176863  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.01659885  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.07418373  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 15] = 0.00241213  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 15] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 15] = 0.019184  / au_2_eV ;

            //Matrix 16:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = -0.03977406  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.00175363  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.01787008  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = -0.03977406  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0065279  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.00576874  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = -0.12214725  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0065279  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.00384038  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = -0.00231661  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.00576874  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = -0.00107321  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = -0.00399243  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.00175363  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.00384038  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = -0.00107321  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.00743378  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.01787008  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = -0.00231661  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = -0.00399243  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.00562151  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 16] = -0.12214725  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 16] = 0.00743378  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 16] = 0.00562151  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 16] = 0.0  / au_2_eV ;

            //Matrix 17:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.03102098  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.00673353  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.03183767  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = -0.01635386  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.02219576  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = -0.014558  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = -0.01074245  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.00673353  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.07088559  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = -0.01635386  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.04530597  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = -0.00999921  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.02219576  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = -0.00999921  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.05646415  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 17] = -0.01074245  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 17] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 17] = -0.02190524  / au_2_eV ;

            //Matrix 18:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = -0.03224536  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = -0.00949155  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.04843608  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.00225872  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = -0.01092102  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.02421777  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.00681468  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = -0.00949155  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.02408195  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.00225872  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = -0.05266955  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.00821411  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = -0.01092102  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.00821411  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0196076  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 18] = 0.00681468  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 18] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 18] = 0.01972858  / au_2_eV ;

            //Matrix 19:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.02228346  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.02228346  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.00286184  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = -0.00175801  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0400158  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.00286184  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = -0.00172272  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = -0.00139876  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = -0.00175801  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.00281512  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = -0.00172272  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = -0.00139876  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.00281512  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.02600811  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 19] = 0.0400158  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 19] = 0.02600811  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 19] = 0.0  / au_2_eV ;

            //Matrix 20:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.00897919  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.01366787  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.04149732  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = -0.00473423  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = -0.01342344  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0840821  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.00154498  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.03523437  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.01366787  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.00154498  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.04258634  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.00122651  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = -0.00473423  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.06368083  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = -0.00534466  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = -0.01342344  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = -0.00534466  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = -0.00993838  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 20] = 0.03523437  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 20] = 0.00122651  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 20] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 20] = 0.06517222  / au_2_eV ;

            //Matrix 21:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = -0.04762007  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.00383849  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = -0.06925937  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = -0.01642414  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.07886478  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.09278602  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.04920817  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.00383849  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.05401473  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = -0.01642414  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = -0.09242341  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0163029  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.07886478  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0163029  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.08807605  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 21] = 0.04920817  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 21] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 21] = 0.05592409  / au_2_eV ;

            //Matrix 22:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.00734702  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.00515335  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = -0.0883017  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.02279118  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.03697325  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = -0.04285892  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = -0.00154002  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.01910082  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.00515335  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = -0.00154002  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.01387786  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.02279118  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = -0.0046263  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = -0.00162272  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.03697325  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = -0.00162272  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = -0.05374124  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 22] = 0.01910082  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 22] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 22] = -0.04857127  / au_2_eV ;

            //Matrix 23:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.01265326  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.00299999  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.00133724  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.16558193  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = -0.00923345  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = -0.03954733  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.07455952  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.00182222  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = -0.00404759  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.00299999  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.00182222  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.00217693  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = -0.00923345  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.05795732  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = -0.00532691  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = -0.03954733  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = -0.00532691  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.07279248  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 23] = 0.00133724  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 23] = -0.00404759  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 23] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 23] = 0.10721237  / au_2_eV ;

            //Matrix 24:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = -0.05891261  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = -0.00446662  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = -0.02258561  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.02977576  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.01904789  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.00897759  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = -0.00446662  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.05020497  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = -0.08073749  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.00978083  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.02977576  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.00978083  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.07175792  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 24] = 0.00897759  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 24] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 24] = 0.02081676  / au_2_eV ;

            //Matrix 25:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.01523759  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.02784554  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.00367381  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.00995579  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.03848751  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = -0.04040859  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.00718377  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.02784554  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.03768853  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.001605  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.00995579  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.08134095  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = -0.01522373  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.03848751  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = -0.01522373  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = -0.05793946  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 25] = 0.00718377  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 25] = 0.001605  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 25] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 25] = -0.05211006  / au_2_eV ;

            //Matrix 26:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.06298634  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.07488392  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.00155457  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.02490458  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.02345776  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.08162395  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = -0.06353845  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = -0.0019606  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = -0.00258349  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.07488392  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = -0.0019606  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.12762941  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.00281735  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.02345776  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.2266666  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = -0.03813577  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.08162395  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = -0.03813577  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = -0.14490265  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 26] = 0.00155457  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 26] = -0.00258349  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 26] = 0.00281735  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 26] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 26] = -0.07986559  / au_2_eV ;

            //Matrix 27:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.01439856  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = -0.00192001  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.01439856  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.00178885  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.00613322  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.00178885  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.00132562  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.00728407  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.00132562  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = -0.00192001  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.00728407  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.00492706  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 27] = 0.00613322  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 27] = 0.00492706  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 27] = 0.0  / au_2_eV ;

            //Matrix 28:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = -0.02693915  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = -0.01079769  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.035783  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = -0.00935391  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.00328403  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = -0.01210925  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = -0.00865102  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = -0.01079769  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.02503434  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = -0.00935391  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = -0.0547027  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.00464929  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.00328403  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.00464929  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.03075676  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 28] = -0.00865102  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 28] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 28] = 0.00163289  / au_2_eV ;

            //Matrix 29:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = -0.0921075  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.04719199  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = -0.06408346  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = -0.02312729  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.01380951  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.05986692  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.00111751  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = -0.01522454  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = -0.04719199  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = -0.00111751  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.07537262  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = -0.0018303  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = -0.02312729  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = -0.21947407  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.04369805  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.01380951  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.04369805  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.25471349  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 29] = -0.01522454  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 29] = -0.0018303  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 29] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 29] = 0.07632584  / au_2_eV ;

            //Matrix 30:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = -0.0017591  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.05651646  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.00195965  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.20694536  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = -0.01468297  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.10085874  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = -0.0017591  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.18001261  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.00267904  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = -0.04603102  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.05651646  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.00267904  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.09048103  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = -0.01468297  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.07564243  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = -0.00706709  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.10085874  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = -0.00706709  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.14680751  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 30] = 0.00195965  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 30] = -0.04603102  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 30] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 30] = 0.17577603  / au_2_eV ;

            //Matrix 31:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.37442957  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = -0.0202115  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = -0.00321373  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.12639056  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = -0.0555621  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = -0.01729615  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = -0.00381004  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = -0.01277169  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = -0.0202115  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = -0.02231451  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.00122069  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = -0.0555621  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.23315431  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = -0.02865154  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = -0.01729615  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = -0.02865154  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.05651686  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 31] = -0.00321373  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 31] = -0.01277169  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 31] = 0.00122069  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 31] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 31] = -0.01387738  / au_2_eV ;

            //Matrix 32:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.09659996  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.00857017  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.06094143  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.05770132  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = -0.0081915  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = -0.05292698  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = -0.00642628  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.00681916  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.00857017  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = -0.00642628  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.42381897  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.00439367  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.05770132  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.13639695  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = -0.02127598  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = -0.0081915  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = -0.02127598  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0849795  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 32] = 0.00681916  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 32] = 0.00439367  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 32] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 32] = -0.0549677  / au_2_eV ;

            //Matrix 33:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = -0.00421776  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = -0.0010213  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = -0.01959208  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = -0.00337991  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = -0.0057256  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = -0.00666662  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 33] = 0.00425324  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = -0.0010213  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = -0.00857158  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = -0.00337991  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = -0.00979591  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = -0.00153129  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = -0.0057256  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = -0.00153129  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = -0.02312994  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 33] = 0.00425324  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 33] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 33] = -0.02775575  / au_2_eV ;

            //Matrix 34:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = -0.00838547  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = -0.00156767  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = 0.00211527  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = -0.00838547  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = -0.00565771  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 34] = 0.01857499  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = -0.00565771  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = -0.00972859  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = -0.00459016  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = -0.00156767  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 34] = 0.00257586  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.00211527  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = -0.00972859  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = -0.00459016  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 34] = 0.03184347  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 34] = 0.01857499  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 34] = 0.00257586  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 34] = 0.03184347  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 34] = 0.0  / au_2_eV ;

            //Matrix 35:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = -0.00625861  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = -0.00761915  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = -0.00124697  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = -0.00336866  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = -0.01523819  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 35] = 0.00449813  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = -0.00136057  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = -0.00124697  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = -0.00367344  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = -0.00146306  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = -0.00336866  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = -0.00146306  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = -0.01591877  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 35] = 0.00449813  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 35] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 35] = -0.02802788  / au_2_eV ;

            //Matrix 36:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = -0.02204115  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = -0.00883842  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = -0.04149733  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = -0.00680998  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.02717995  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 36] = 0.06731384  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = -0.00883842  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = -0.00435389  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 36] = 0.00108336  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = -0.04013906  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.00445954  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = -0.00680998  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.00445954  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = -0.00571208  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 36] = 0.06731384  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 36] = 0.00108336  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 36] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 36] = 0.03907974  / au_2_eV ;

            //Matrix 37:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.02068061  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.00442597  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.01115659  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = -0.0016119  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.00722506  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = -0.00843525  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 37] = 0.00512049  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.00442597  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0221773  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = -0.0016119  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = -0.02557868  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.00722506  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 37] = 0.00512049  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 37] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 37] = -0.00353776  / au_2_eV ;

            //Matrix 38:
            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = -0.00952405  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.00216377  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.00572697  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[0 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 38] = 0.00152505  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.0200004  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = -0.00874246  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.00941607  / au_2_eV ;

            setm->c_lvcm[1 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.00216377  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.17476473  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[2 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 38] = 0.09498998  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.00572697  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = -0.00802736  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[3 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = -0.00874246  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = -0.00340227  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.00340416  / au_2_eV ;

            setm->c_lvcm[4 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.00941607  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.00340416  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.01728002  / au_2_eV ;

            setm->c_lvcm[5 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 0 * setm->N_mode_lvcm + 38] = 0.00152505  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 1 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 2 * setm->N_mode_lvcm + 38] = 0.09498998  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 3 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 4 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 5 * setm->N_mode_lvcm + 38] = 0.0  / au_2_eV ;

            setm->c_lvcm[6 * setm->Nstate_lvcm * setm->N_mode_lvcm + 6 * setm->N_mode_lvcm + 38] = 0.09626052  / au_2_eV ;


            break;
     }

     
}

void sample_dnalvcm(double *P, double *R, struct set_host *setm) {
    int j;
    double x2;

    sample_LVCM(P, R, setm);
    

}

void V_dnalvcm(double *R, double complex *H, int forcetype, struct set_host *setm) {
    
    V_LVCM(R, H, forcetype, setm);
    
}

void dV_dnalvcm(double *R, double complex *dH, int forcetype, struct set_host *setm) {
    
    dV_LVCM(R, dH, forcetype, setm);

}

void nucforce_dnalvcm(double *R, double *nf, struct set_host *setm) {

    nucforce_LVCM(R, nf, setm);

}



