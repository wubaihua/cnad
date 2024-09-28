#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"
#include <string.h>
#include <stdio.h>
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif




void parameter_SEM(double *mass, struct set_host *setm) ;

void sample_SEM(double *P, double *R, double beta, struct set_host *setm) ;

void V_SEM(double *R, double *H, int forcetype, struct set_host *setm) ;

void dV_SEM(double *R, double *dH, int forcetype, struct set_host *setm) ;

void nucforce_SEM(double *R, double *nf, struct set_host *setm);