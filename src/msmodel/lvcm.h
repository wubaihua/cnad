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




void parameter_LVCM(double *mass, struct set_host *setm);

void sample_LVCM(double *P, double *R, struct set_host *setm);

void V_LVCM(double *R, double complex *H, int forcetype, struct set_host *setm) ;

void dV_LVCM(double *R, double complex *dH, int forcetype, struct set_host *setm) ;

void nucforce_LVCM(double *R, double *nf, struct set_host *setm);

