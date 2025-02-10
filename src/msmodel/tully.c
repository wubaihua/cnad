#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif


// int type_tully;
// double A_tully, B_tully, C_tully, D_tully, E_tully, A2_tully, f_tully, Z_tully;
// double R0_tully, P0_tully, gamma_tully;

void parameter_tully(double *mass, struct set_host *setm) {
    
    
    mass[0] = 2000;

    switch (setm->type_tully) {
        case 1:
            setm->A_tully = 0.01;
            setm->B_tully = 1.6;
            setm->C_tully = 0.005;
            setm->D_tully = 1;
            setm->R0_tully = -3.8;
            break;
        case 2:
            setm->A_tully = 0.1;
            setm->B_tully = 0.28;
            setm->C_tully = 0.015;
            setm->D_tully = 0.06;
            setm->E_tully = 0.05;
            setm->R0_tully = -10;
            break;
        case 3:
            setm->A_tully = 6E-4;
            setm->B_tully = 0.1;
            setm->C_tully = 0.9;
            setm->R0_tully = -13;
            break;
        case 4:
            setm->A_tully = 0.04;
            setm->A2_tully = 0.01;
            setm->B_tully = 1.0;
            setm->C_tully = 0.005;
            setm->D_tully = 1;
            setm->f_tully = 0.7;
            mass[0] = 1980;
            setm->R0_tully = -5;
            break;
        case 5:
            setm->A_tully = 6E-4;
            setm->B_tully = 0.1;
            setm->C_tully = 0.9;
            setm->Z_tully = 4;
            setm->R0_tully = -27.5;
            break;
        case 6:
            mass[0] = 2000.0;
            setm->R0_tully = 0.0;
            break;

    }
}

void sample_tully(double *P, double *R, struct set_host *setm) {
    double x2;

    if(setm->if_flighttime_tully == 0){
        switch (setm->type_tully) {
        case 1:
        case 2:
        case 3:
            setm->gamma_tully = 1;
            break;
        case 4:
            setm->gamma_tully = 0.25;
            break;
        case 5:
            setm->gamma_tully = 1.0 / 8;
            break;
        case 6:

            break;
        }
    }
    

    box_muller(&P[0], &x2, sqrt(setm->gamma_tully / 2.0), setm->P0_tully);
    box_muller(&R[0], &x2, sqrt(1.0 / (setm->gamma_tully * 2.0)), setm->R0_tully);


    if(setm->if_flighttime_tully == 1){
        R[0] = -1.0 * setm->Xb_tully;
    }

    // P[0] = setm->P0_tully;
}

void V_tully(double *R, double *H, struct set_host *setm) {
    switch (setm->type_tully) {
        case 1:
            if (R[0] > 0) {
                H[0] = setm->A_tully * (1 - exp(-setm->B_tully * R[0]));
            } else {
                H[0] = -setm->A_tully * (1 - exp(setm->B_tully * R[0]));
            }
            H[3] = -H[0];
            H[1] = setm->C_tully * exp(-setm->D_tully * R[0] * R[0]);
            H[2] = H[1];
            break;
        case 2:
            H[0] = 0;
            H[3] = -setm->A_tully * exp(-setm->B_tully * R[0] * R[0]) + setm->E_tully;
            H[1] = setm->C_tully * exp(-setm->D_tully * R[0] * R[0]);
            H[2] = H[1];
            break;
        case 3:
            H[0] = -setm->A_tully;
            H[3] = setm->A_tully;
            if (R[0] > 0) {
                H[1] = setm->B_tully * (1 - (exp(-setm->C_tully * R[0]) - 1));
            } else {
                H[1] = setm->B_tully * (1 + (exp(setm->C_tully * R[0]) - 1));
            }
            H[2] = H[1];
            break;
        case 4:
            H[0] = setm->A_tully * (1 + tanh(setm->B_tully * R[0]));
            H[3] = setm->A2_tully * (1 - tanh(setm->B_tully * R[0]));
            H[1] = setm->C_tully * exp(-setm->D_tully * (R[0] + setm->f_tully) * (R[0] + setm->f_tully));
            H[2] = H[1];
            break;
        case 5:
            H[0] = setm->A_tully;
            H[3] = -setm->A_tully;
            if (R[0] < -setm->Z_tully) {
                H[1] = -setm->B_tully * exp(setm->C_tully * (R[0] - setm->Z_tully)) + setm->B_tully * exp(setm->C_tully * (R[0] + setm->Z_tully));
            } else if (R[0] > setm->Z_tully) {
                H[1] = setm->B_tully * exp(-setm->C_tully * (R[0] - setm->Z_tully)) - setm->B_tully * exp(-setm->C_tully * (R[0] + setm->Z_tully));
            } else {
                H[1] = -setm->B_tully * exp(setm->C_tully * (R[0] - setm->Z_tully)) - setm->B_tully * exp(-setm->C_tully * (R[0] + setm->Z_tully)) + 2 * setm->B_tully;
            }
            H[2] = H[1];
            break;
        case 6:
            if (R[0] < 2 * M_PI){
                H[0] = 0.5 * setm->A_tully;
                H[1] = 0.0;
                H[2] = 0.0;
                H[3] = -0.5 * setm->A_tully;
            } else {
                H[0] = 0.5 * setm->A_tully * cos(setm->B_tully * R[0]);
                H[1] = 0.5 * setm->A_tully * sin(setm->B_tully * R[0]);
                H[2] = 0.5 * setm->A_tully * sin(setm->B_tully * R[0]);
                H[3] = -0.5 * setm->A_tully * cos(setm->B_tully * R[0]);
            }
            break;
    }
}

void dV_tully(double *R, double *dH, struct set_host *setm) {
    switch (setm->type_tully) {
        case 1:
            if (R[0] > 0) {
                dH[0] = setm->A_tully * setm->B_tully * exp(-setm->B_tully * R[0]);
            } else {
                dH[0] = setm->A_tully * setm->B_tully * exp(setm->B_tully * R[0]);
            }
            dH[3] = -dH[0];
            dH[1] = -setm->C_tully * 2 * setm->D_tully * R[0] * exp(-setm->D_tully * R[0] * R[0]);
            dH[2] = dH[1];
            break;
        case 2:
            dH[0] = 0;
            dH[3] = 2 * setm->A_tully * setm->B_tully * R[0] * exp(-setm->B_tully * R[0] * R[0]);
            dH[1] = -setm->C_tully * 2 * setm->D_tully * R[0] * exp(-setm->D_tully * R[0] * R[0]);
            dH[2] = dH[1];
            break;
        case 3:
            dH[0] = 0;
            dH[3] = 0;
            if (R[0] > 0) {
                dH[1] = setm->B_tully * setm->C_tully * exp(-setm->C_tully * R[0]);
            } else {
                dH[1] = setm->B_tully * setm->C_tully * exp(setm->C_tully * R[0]);
            }
            dH[2] = dH[1];
            break;
        case 4:
            dH[0] = setm->A_tully * setm->B_tully / cosh(setm->B_tully * R[0]) / cosh(setm->B_tully * R[0]);
            dH[3] = -setm->A2_tully * setm->B_tully / cosh(setm->B_tully * R[0]) / cosh(setm->B_tully * R[0]);
            dH[1] = -setm->C_tully * 2 * setm->D_tully * (R[0] + setm->f_tully) * exp(-setm->D_tully * (R[0] + setm->f_tully) * (R[0] + setm->f_tully));
            dH[2] = dH[1];
            break;
        case 5:
            dH[0] = 0;
            dH[3] = 0;
            if (R[0] < -setm->Z_tully) {
                dH[1] = -setm->B_tully * exp(setm->C_tully * (R[0] - setm->Z_tully)) * setm->C_tully + setm->B_tully * exp(setm->C_tully * (R[0] + setm->Z_tully)) * setm->C_tully;
            } else if (R[0] > setm->Z_tully) {
                dH[1] = -setm->C_tully * setm->B_tully * exp(-setm->C_tully * (R[0] - setm->Z_tully)) + setm->C_tully * setm->B_tully * exp(-setm->C_tully * (R[0] + setm->Z_tully));
            } else {
                dH[1] = -setm->B_tully * exp(setm->C_tully * (R[0] - setm->Z_tully)) * setm->C_tully + setm->B_tully * exp(-setm->C_tully * (R[0] + setm->Z_tully)) * setm->C_tully;
            }
            dH[2] = dH[1];
            break;
        case 6:
            if (R[0] < 2 * M_PI){
                dH[0] = 0.0;
                dH[1] = 0.0;
                dH[2] = 0.0;
                dH[3] = 0.0;
            } else {
                dH[0] = -0.5 * setm->A_tully * setm->B_tully * sin(setm->B_tully * R[0]);
                dH[1] = 0.5 * setm->A_tully * setm->B_tully * cos(setm->B_tully * R[0]);
                dH[2] = 0.5 * setm->A_tully * setm->B_tully * cos(setm->B_tully * R[0]);
                dH[3] = 0.5 * setm->A_tully * setm->B_tully * sin(setm->B_tully * R[0]);
            }
            break;
    }
}