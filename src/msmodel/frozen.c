#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif



void parameter_frozen(double *mass, struct set_host *setm) {
    
    
    mass[0] = 1;

    
}

void sample_frozen(double *P, double *R, struct set_host *setm) {
   
    P[0]=0;
    R[0]=0;

}

void V_frozen(double *R, double *H, struct set_host *setm) {
    memcpy(H,setm->Hele_frozen,setm->Nstate_frozen * setm->Nstate_frozen * sizeof(double));
}

void dV_frozen(double *R, double *dH, struct set_host *setm) {
    memset(dH,0,setm->Nstate_frozen * setm->Nstate_frozen * sizeof(double));
}




void cfweight_frozen(double *w0, double *wt, struct set_host *setm) {
    int i, j;

    // for (i = 0; i < n; i++) {
    //     w0[i] = 0.0;
    //     wt[i] = 0.0;
    // }

    // for (i = 0; i < setm->Nstate_SEMdp - 1; i++) {
    //     w0[i * setm->Nstate_SEMdp + (setm->Nstate_SEMdp - 1)] = setm->dipole_SEMdp[i];
    //     wt[i * setm->Nstate_SEMdp + (setm->Nstate_SEMdp - 1)] = setm->dipole_SEMdp[i];
    //     wt[(setm->Nstate_SEMdp - 1) * setm->Nstate_SEMdp + i] = setm->dipole_SEMdp[i];
    // }


    
    // double a1[] = {1.0,2.0,4.0,
    //                3.0,2.0,8.0,
    //                6.0,7.0,4.0};
    // double a2[] = {0.0,1.0,3.0,
    //                -2.0,0.0,1.0,
    //                9.0,-7.0,0.0};

    // double a1[] = {1.0,0.0,0.0,
    //                0.0,2.0,0.0,
    //                0.0,0.0,4.0};
    // double a2[] = {0.0,1.0,3.0,
    //                -2.0,0.0,1.0,
    //                9.0,-7.0,0.0};

     double a1[] = {1.0,  2.0, 4.0, 6.0, 7.0,
                    3.0,  2.0, 8.0, 7.0, 8.0,
                    6.0,  7.0, 4.0, 5.0, 6.0,
                    6.0, -7.0, 4.0, 5.0, 3.0,
                    -7.0, -1.0, 2.0, 9.0, 7.0};
                                  
    double a2[] = {0.0,  1.0, 3.0,  6.0, -7.0,\
                    -2.0,  0.0, 1.0, -7.0,  3.0,\
                    9.0, -7.0, 0.0,  5.0,  2.0,\
                    6.0,  7.0, 4.0,  0.0,  5.0,\
                    7.0,  1.0, 6.0, -7.0,  0.0};             


    for (i = 0; i < setm->Nstate_frozen*setm->Nstate_frozen; i++) {
        w0[i] = a1[i];
        wt[i] = a2[i];
    }
}
