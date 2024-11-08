#include <stdio.h>
#include <math.h>
#include "gmath.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif



void parameter_frozen(double *mass, struct set_host *setm);

void sample_frozen(double *P, double *R, struct set_host *setm);

void V_frozen(double *R, double *H, struct set_host *setm);

void dV_frozen(double *R, double *dH, struct set_host *setm);