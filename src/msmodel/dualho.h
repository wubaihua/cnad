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

void parameter_dualho(double *mass, struct set_host *setm);

void sample_dualho(double *P, double *R, struct set_host *setm) ;

void V_dualho(double *R, double *H, struct set_host *setm);
void dV_dualho(double *R, double *dH, struct set_host *setm);