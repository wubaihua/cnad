#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"  // Assuming you have a math.h header for functions like box_muller
// #include "constant.h"
#include <string.h>
#include <stdio.h>
#include "cJSON.h"

extern int forcetype;
extern char msmodelname[200];


// Spin-Boson Model parameters
extern int N_bath_SBM, bathtype; // bathtype=1 for Ohmic; bathtype=2 for Debye
extern double eps_SBM, delta_SBM, alpha_SBM, omega_c_SBM, lambda_SBM, s_SBM;
void readinp_SBM(cJSON *item, int *Ndof1, int *Ndof2, int *Nstate);








void readinp_msmodel(cJSON *json, int *Ndof1, int *Ndof2, int *Nstate);

