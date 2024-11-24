#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif


// int type_dualho;
// double A_dualho, B_dualho, C_dualho, D_dualho, E_dualho, A2_dualho, f_dualho, Z_dualho;
// double R0_dualho, P0_dualho, gamma_dualho;

void parameter_dualho(double *mass, struct set_host *setm) {
    
    
    mass[0] = 20000;

    // switch (setm->type_dualho) {
    //     case 1:
    //         setm->A_dualho = 0.01;
    //         setm->B_dualho = 1.6;
    //         setm->C_dualho = 0.005;
    //         setm->D_dualho = 1;
    //         setm->R0_dualho = -3.8;
    //         break;
    //     case 2:
    //         setm->A_dualho = 0.1;
    //         setm->B_dualho = 0.28;
    //         setm->C_dualho = 0.015;
    //         setm->D_dualho = 0.06;
    //         setm->E_dualho = 0.05;
    //         setm->R0_dualho = -10;
    //         break;
    //     case 3:
    //         setm->A_dualho = 6E-4;
    //         setm->B_dualho = 0.1;
    //         setm->C_dualho = 0.9;
    //         setm->R0_dualho = -13;
    //         break;
    //     case 4:
    //         setm->A_dualho = 0.04;
    //         setm->A2_dualho = 0.01;
    //         setm->B_dualho = 1.0;
    //         setm->C_dualho = 0.005;
    //         setm->D_dualho = 1;
    //         setm->f_dualho = 0.7;
    //         mass[0] = 1980;
    //         setm->R0_dualho = -5;
    //         break;
    //     case 5:
    //         setm->A_dualho = 6E-4;
    //         setm->B_dualho = 0.1;
    //         setm->C_dualho = 0.9;
    //         setm->Z_dualho = 4;
    //         setm->R0_dualho = -27.5;
    //         break;
    // }
}

void sample_dualho(double *P, double *R, struct set_host *setm) {
    double x2;

    switch (setm->type_dualho) {
        case 1:
            box_muller(&P[0], &x2, 1.0/(2*0.2236), 0.0);
            box_muller(&R[0], &x2, 0.2236, 2.0);
    }

    
}

void V_dualho(double *R, double *H, struct set_host *setm) {
    switch (setm->type_dualho) {
        case 1:
            
            H[0] = 0.01 * (R[0] - 6.0) * (R[0] - 6.0);
            H[3] = 0.01 * (R[0] - 2.0) * (R[0] - 2.0) + 0.01;
            H[1] = 0.01 * exp(-3 * (R[0] - 3.875) * (R[0] - 3.875));
            H[2] = H[1];
            break;
        // case 2:
        //     H[0] = 0;
        //     H[3] = -setm->A_dualho * exp(-setm->B_dualho * R[0] * R[0]) + setm->E_dualho;
        //     H[1] = setm->C_dualho * exp(-setm->D_dualho * R[0] * R[0]);
        //     H[2] = H[1];
        //     break;
        // case 3:
        //     H[0] = -setm->A_dualho;
        //     H[3] = setm->A_dualho;
        //     if (R[0] > 0) {
        //         H[1] = setm->B_dualho * (1 - (exp(-setm->C_dualho * R[0]) - 1));
        //     } else {
        //         H[1] = setm->B_dualho * (1 + (exp(setm->C_dualho * R[0]) - 1));
        //     }
        //     H[2] = H[1];
        //     break;
        // case 4:
        //     H[0] = setm->A_dualho * (1 + tanh(setm->B_dualho * R[0]));
        //     H[3] = setm->A2_dualho * (1 - tanh(setm->B_dualho * R[0]));
        //     H[1] = setm->C_dualho * exp(-setm->D_dualho * (R[0] + setm->f_dualho) * (R[0] + setm->f_dualho));
        //     H[2] = H[1];
        //     break;
        // case 5:
        //     H[0] = setm->A_dualho;
        //     H[3] = -setm->A_dualho;
        //     if (R[0] < -setm->Z_dualho) {
        //         H[1] = -setm->B_dualho * exp(setm->C_dualho * (R[0] - setm->Z_dualho)) + setm->B_dualho * exp(setm->C_dualho * (R[0] + setm->Z_dualho));
        //     } else if (R[0] > setm->Z_dualho) {
        //         H[1] = setm->B_dualho * exp(-setm->C_dualho * (R[0] - setm->Z_dualho)) - setm->B_dualho * exp(-setm->C_dualho * (R[0] + setm->Z_dualho));
        //     } else {
        //         H[1] = -setm->B_dualho * exp(setm->C_dualho * (R[0] - setm->Z_dualho)) - setm->B_dualho * exp(-setm->C_dualho * (R[0] + setm->Z_dualho)) + 2 * setm->B_dualho;
        //     }
        //     H[2] = H[1];
        //     break;
    }
}

void dV_dualho(double *R, double *dH, struct set_host *setm) {
    switch (setm->type_dualho) {
        case 1:
            
            dH[0] = 0.01 * 2 * (R[0] - 6.0);
            dH[3] = 0.01 * 2 * (R[0] - 2.0);
            dH[1] = -0.01 * exp(-3 * (R[0] - 3.875) * (R[0] - 3.875)) * 6 * (R[0] - 3.875);
            dH[2] = dH[1];
            break;
        // case 2:
        //     dH[0] = 0;
        //     dH[3] = 2 * setm->A_dualho * setm->B_dualho * R[0] * exp(-setm->B_dualho * R[0] * R[0]);
        //     dH[1] = -setm->C_dualho * 2 * setm->D_dualho * R[0] * exp(-setm->D_dualho * R[0] * R[0]);
        //     dH[2] = dH[1];
        //     break;
        // case 3:
        //     dH[0] = 0;
        //     dH[3] = 0;
        //     if (R[0] > 0) {
        //         dH[1] = setm->B_dualho * setm->C_dualho * exp(-setm->C_dualho * R[0]);
        //     } else {
        //         dH[1] = setm->B_dualho * setm->C_dualho * exp(setm->C_dualho * R[0]);
        //     }
        //     dH[2] = dH[1];
        //     break;
        // case 4:
        //     dH[0] = setm->A_dualho * setm->B_dualho / cosh(setm->B_dualho * R[0]) / cosh(setm->B_dualho * R[0]);
        //     dH[3] = -setm->A2_dualho * setm->B_dualho / cosh(setm->B_dualho * R[0]) / cosh(setm->B_dualho * R[0]);
        //     dH[1] = -setm->C_dualho * 2 * setm->D_dualho * (R[0] + setm->f_dualho) * exp(-setm->D_dualho * (R[0] + setm->f_dualho) * (R[0] + setm->f_dualho));
        //     dH[2] = dH[1];
        //     break;
        // case 5:
        //     dH[0] = 0;
        //     dH[3] = 0;
        //     if (R[0] < -setm->Z_dualho) {
        //         dH[1] = -setm->B_dualho * exp(setm->C_dualho * (R[0] - setm->Z_dualho)) * setm->C_dualho + setm->B_dualho * exp(setm->C_dualho * (R[0] + setm->Z_dualho)) * setm->C_dualho;
        //     } else if (R[0] > setm->Z_dualho) {
        //         dH[1] = -setm->C_dualho * setm->B_dualho * exp(-setm->C_dualho * (R[0] - setm->Z_dualho)) + setm->C_dualho * setm->B_dualho * exp(-setm->C_dualho * (R[0] + setm->Z_dualho));
        //     } else {
        //         dH[1] = -setm->B_dualho * exp(setm->C_dualho * (R[0] - setm->Z_dualho)) * setm->C_dualho + setm->B_dualho * exp(-setm->C_dualho * (R[0] + setm->Z_dualho)) * setm->C_dualho;
        //     }
        //     dH[2] = dH[1];
        //     break;
    }
}