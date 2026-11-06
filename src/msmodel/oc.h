#ifndef OC_H
#define OC_H

#include <complex.h>
#include "def_host.h"

void parameter_OC(double *mass, struct set_host *setm);
void sample_OC(double *P, double *R, double beta, struct set_host *setm);
void V_OC(double *R, double *H, int forcetype, struct set_host *setm);
void dV_OC(double *R, double *dH, int forcetype, struct set_host *setm);
void nucforce_OC(double *R, double *nf, struct set_host *setm);
void cfweight_OC(double *w0, double *wt, double beta, double *R, double *P,
                 double *Rt, struct set_host *setm);
void mqcequden_OC(double *rho, double beta, double *R, double *P,
                  struct set_host *setm);
/* U and E must already be phase-matched by the shared cal_NACV setup. */
void cal_NACV_oc(const double complex *U, const double complex *dV,
                 const double *E, double complex *dv_adia, double complex *nac,
                 struct set_host *setm);

#endif
