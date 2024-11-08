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