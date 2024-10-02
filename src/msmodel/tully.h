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

void parameter_tully(double *mass, struct set_host *setm);

void sample_tully(double *P, double *R, struct set_host *setm);

void V_tully(double *R, double *H, struct set_host *setm);

void dV_tully(double *R, double *dH, struct set_host *setm);