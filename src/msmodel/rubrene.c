#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"  // Assuming you have a math.h header for functions like box_muller
// #include "constant.h"
#include <string.h>
#include <stdio.h>
// #include "cJSON.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif

// Constants and parameters
// int setm->Nstate_rubrene, setm->N_mode_rubrene, setm->Nmole_rubrene;
// double setm->Vc_rubrene, setm->factor_freq_rubrene;
// double *setm->omega_rubrene, *setm->lambda_rubrene, *setm->g_rubrene, *setm->mass_rubrene;

// void readinp_rubrene(FILE *idinp, int *Ndof1, int *Ndof2, int *Nstate) {
//     char c50[50];
    
//     fscanf(idinp, "%s", c50); // Read character string
//     fscanf(idinp, "%d", &setm->Nmole_rubrene);
//     fscanf(idinp, "%s", c50);
//     fscanf(idinp, "%lf", &setm->Vc_rubrene);
//     fscanf(idinp, "%s", c50);
//     fscanf(idinp, "%lf", &setm->factor_freq_rubrene);

//     setm->Vc_rubrene = setm->Vc_rubrene / au_2_eV;

//     setm->Nstate_rubrene = setm->Nmole_rubrene;
//     setm->N_mode_rubrene = 9;

//     *Ndof1 = setm->Nmole_rubrene;
//     *Ndof2 = setm->N_mode_rubrene;
//     *Nstate = setm->Nstate_rubrene;
// }




void parameter_rubrene(double *mass, struct set_host *setm) {
    int i;
    
    
    for (i = 0; i < setm->Nmole_rubrene * setm->N_mode_rubrene; i++) {
        mass[i] = 1.0 / (setm->factor_freq_rubrene * setm->factor_freq_rubrene);
    }

    setm->omega_rubrene = (double *)malloc(setm->N_mode_rubrene * sizeof(double));
    setm->lambda_rubrene = (double *)malloc(setm->N_mode_rubrene * sizeof(double));
    setm->g_rubrene = (double *)malloc(setm->N_mode_rubrene * sizeof(double));
    setm->mass_rubrene = (double *)malloc(setm->Nmole_rubrene * setm->N_mode_rubrene * sizeof(double));

    double omega_values[] = {10.0, 27.0, 78.0, 124.0, 149.0, 167.0, 169.0, 190.0, 198.0};
    double g_values[] = {0.96, 0.38, 0.25, 0.20, 0.15, 0.31, 0.13, 0.20, 0.31};

    for (i = 0; i < setm->N_mode_rubrene; i++) {
        setm->omega_rubrene[i] = omega_values[i] * 0.001 / au_2_eV * setm->factor_freq_rubrene;
        setm->g_rubrene[i] = g_values[i] / sqrt(setm->factor_freq_rubrene);
        setm->lambda_rubrene[i] = setm->g_rubrene[i] * setm->g_rubrene[i] * setm->omega_rubrene[i];
    }

    for (i = 0; i < setm->Nmole_rubrene * setm->N_mode_rubrene; i++) {
        setm->mass_rubrene[i] = mass[i];
    }
}

void sample_rubrene(double *P, double *R, double beta, struct set_host *setm) {
    int i, j;
    double x2;

    

    for (i = 0; i < setm->Nmole_rubrene; i++) {
        for (j = 0; j < setm->N_mode_rubrene; j++) {
            if (beta > 99999) {
                if(setm->if_classical == 0){
                    box_muller(&P[i * setm->N_mode_rubrene + j], &x2, sqrt(0.5 * hbar * setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j]), 0.0);
                    box_muller(&R[i * setm->N_mode_rubrene + j], &x2, sqrt(0.5 * hbar / (setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j])), 0.0);
                } else {
                    P[i * setm->N_mode_rubrene + j] = 0;
                    R[i * setm->N_mode_rubrene + j] = 0;
                }
                
            } else {
                if(setm->if_classical == 0){
                    box_muller(&P[i * setm->N_mode_rubrene + j], &x2, sqrt(0.5 * hbar * setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j] / tanh(0.5 * beta * hbar * setm->omega_rubrene[j])), 0.0);
                    box_muller(&R[i * setm->N_mode_rubrene + j], &x2, sqrt(0.5 * hbar / (tanh(0.5 * beta * hbar * setm->omega_rubrene[j]) * setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j])), 0.0);
                } else {
                    box_muller(&P[i * setm->N_mode_rubrene + j], &x2, sqrt(setm->mass_rubrene[i * setm->N_mode_rubrene + j]  / beta), 0.0);
                    box_muller(&R[i * setm->N_mode_rubrene + j], &x2, sqrt(1.0 / (beta * setm->omega_rubrene[j] * setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j])), 0.0);
                }
            }
        }
    }
}

void V_rubrene(double *R, double *H, int forcetype, struct set_host *setm) {
    int i, j, k;

    
    for (i = 0; i < setm->Nstate_rubrene * setm->Nstate_rubrene; i++) {
        H[i] = 0.0;
    }

    for (i = 0; i < setm->Nstate_rubrene; i++) {
        double sum_val = 0.0;
        for (j = 0; j < setm->N_mode_rubrene; j++) {
            sum_val += sqrt(2 * setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j]) * setm->g_rubrene[j] * setm->omega_rubrene[j] * R[i * setm->N_mode_rubrene + j];
        }
        H[i * setm->Nstate_rubrene + i] = sum_val;

        if (i < setm->Nstate_rubrene - 1) {
            H[i * setm->Nstate_rubrene + (i + 1)] = setm->Vc_rubrene;
            H[(i + 1) * setm->Nstate_rubrene + i] = setm->Vc_rubrene;
        } else {
            H[0 * setm->Nstate_rubrene + (setm->Nstate_rubrene - 1)] = setm->Vc_rubrene;
            H[(setm->Nstate_rubrene - 1) * setm->Nstate_rubrene + 0] = setm->Vc_rubrene;
        }
    }

    if (forcetype == 0) {
        double Vnuc = 0;
        for (i = 0; i < setm->Nmole_rubrene; i++) {
            for (j = 0; j < setm->N_mode_rubrene; j++) {
                Vnuc += 0.5 * setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j] * setm->omega_rubrene[j] * R[i * setm->N_mode_rubrene + j] * R[i * setm->N_mode_rubrene + j];
            }
        }

        for (i = 0; i < setm->Nstate_rubrene; i++) {
            H[i * setm->Nstate_rubrene + i] += Vnuc;
        }
    }
}

void dV_rubrene(double *R, double *dH, int forcetype, struct set_host *setm) {
    int i, j;

    int dH_size = setm->Nstate_rubrene * setm->Nstate_rubrene * setm->Nmole_rubrene * setm->N_mode_rubrene;
   
    for (i = 0; i < dH_size; i++) {
        dH[i] = 0.0;
    }

    for (i = 0; i < setm->Nstate_rubrene; i++) {
        for (j = 0; j < setm->N_mode_rubrene; j++) {
            dH[i * setm->Nstate_rubrene * setm->Nmole_rubrene * setm->N_mode_rubrene 
                + i * setm->Nmole_rubrene * setm->N_mode_rubrene 
                + i * setm->N_mode_rubrene + j] = sqrt(2 * setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j]) * setm->g_rubrene[j] * setm->omega_rubrene[j];
        }
    }

    if (forcetype == 0) {
        double *force = (double *)malloc(setm->Nmole_rubrene * setm->N_mode_rubrene * sizeof(double));
        for (int i = 0; i < setm->Nmole_rubrene; i++) {
            for (int j = 0; j < setm->N_mode_rubrene; j++) {
                force[i * setm->N_mode_rubrene + j] = setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j] * setm->omega_rubrene[j] * R[i * setm->N_mode_rubrene + j];
            }
        }

        for (int i = 0; i < setm->Nmole_rubrene; i++) {
            for (int j = 0; j < setm->Nmole_rubrene; j++) {
                for (int k = 0; k < setm->N_mode_rubrene; k++) {
                    dH[  i * setm->Nmole_rubrene * setm->Nmole_rubrene * setm->N_mode_rubrene 
                    + i * setm->Nmole_rubrene * setm->N_mode_rubrene + j * setm->N_mode_rubrene + k] += force[j * setm->N_mode_rubrene + k];
                }
            }
        }
        free(force); 
    }
}

void nucforce_rubrene(double *R, double *nf, struct set_host *setm) {
    int i, j;

    
    for (i = 0; i < setm->Nmole_rubrene * setm->N_mode_rubrene; i++) {
        nf[i] = 0.0;
    }

    for (i = 0; i < setm->Nmole_rubrene; i++) {
        for (j = 0; j < setm->N_mode_rubrene; j++) {
            nf[i * setm->N_mode_rubrene + j] = setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j] * setm->omega_rubrene[j] * R[i * setm->N_mode_rubrene + j];
        }
    }
}




void cfweight_rubrene(double *w0, double *wt, double beta, double *R, double *P, struct set_host *setm) {
    int i, j;
    double *rho, *Heff, *E, *C, *expe;
    double *tempdm1, *tempdm2; 


    rho = (double *)malloc(setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
    Heff = (double *)malloc(setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
    E = (double *)malloc(setm->Nstate_rubrene * sizeof(double));
    C = (double *)malloc(setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
    expe = (double *)malloc(setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
    if(setm->if_classical == 0){
        
        
        memset(Heff, 0, setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
 
        for (i = 0; i < setm->Nstate_rubrene; i++) {
            if (i < setm->Nstate_rubrene - 1) {
                Heff[i * setm->Nstate_rubrene + (i + 1)] = setm->Vc_rubrene;
                Heff[(i + 1) * setm->Nstate_rubrene + i] = setm->Vc_rubrene;
            } else {
                Heff[0 * setm->Nstate_rubrene + (setm->Nstate_rubrene - 1)] = setm->Vc_rubrene;
                Heff[(setm->Nstate_rubrene - 1) * setm->Nstate_rubrene + 0] = setm->Vc_rubrene;
            }
        }
 
        double sum1 = 0.0;
        for (i = 0; i < setm->N_mode_rubrene; i++){
            sum1 += setm->lambda_rubrene[i];
        }
        for (i = 0; i < setm->Nstate_rubrene * setm->Nstate_rubrene; i++) {
            Heff[i] *= exp(- sum1 * beta / 3.0);
        }
 
        dia_symmat(setm->Nstate_rubrene, Heff, E, C);
 
        memset(expe, 0, setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        for (i = 0; i < setm->Nstate_rubrene; i++) {
            expe[i * setm->Nstate_rubrene + i] = exp(-beta * E[i]);
        }
 
       
 
        tempdm1 = (double *)malloc(setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        tempdm2 = (double *)malloc(setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        transpose(C,tempdm1,setm->Nstate_rubrene);
        dd_matmul(C,expe,tempdm2,setm->Nstate_rubrene,setm->Nstate_rubrene,setm->Nstate_rubrene);
        dd_matmul(tempdm2,tempdm1,rho,setm->Nstate_rubrene,setm->Nstate_rubrene,setm->Nstate_rubrene);
 
        double sum = 0.0;
        for (i = 0; i < setm->Nstate_rubrene; i++) {
            sum += rho[i * setm->Nstate_rubrene + i];
        }
 
        
        for (i = 0; i < setm->Nstate_rubrene * setm->Nstate_rubrene; i++) {
            rho[i] /= sum;
        }
 
        memset(w0, 0, setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        memset(wt, 0, setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
 
        for (i = 0; i < setm->Nstate_rubrene; i++) {
            if (i == 0) {
                wt[(i + 1) * setm->Nstate_rubrene + i] = 1;
                wt[i * setm->Nstate_rubrene + (i + 1)] = -1;
                for (j = 0; j < setm->Nstate_rubrene; j++) {
                    w0[i * setm->Nstate_rubrene + j] = rho[(setm->Nstate_rubrene - 1) * setm->Nstate_rubrene + j] - rho[(i + 1) * setm->Nstate_rubrene + j];
                }
            } else if (i == setm->Nstate_rubrene - 1) {
                wt[0 * setm->Nstate_rubrene + i] = 1;
                wt[i * setm->Nstate_rubrene + 0] = -1;
                for (j = 0; j < setm->Nstate_rubrene; j++) {
                    w0[i * setm->Nstate_rubrene + j] = rho[(i - 1) * setm->Nstate_rubrene + j] - rho[0 * setm->Nstate_rubrene + j];
                }
            } else {
                wt[(i + 1) * setm->Nstate_rubrene + i] = 1;
                wt[i * setm->Nstate_rubrene + (i + 1)] = -1;
                for (j = 0; j < setm->Nstate_rubrene; j++) {
                    w0[i * setm->Nstate_rubrene + j] = rho[(i - 1) * setm->Nstate_rubrene + j] - rho[(i + 1) * setm->Nstate_rubrene + j];
                }
            }
        }
 
        

    } else {
        memset(Heff, 0, setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
 
        for (i = 0; i < setm->Nstate_rubrene; i++) {
            double sum_val = 0.0;
            for (j = 0; j < setm->N_mode_rubrene; j++) {
                sum_val += sqrt(2 * setm->mass_rubrene[i * setm->N_mode_rubrene + j] * setm->omega_rubrene[j]) * setm->g_rubrene[j] * setm->omega_rubrene[j] * R[i * setm->N_mode_rubrene + j];
            }
            Heff[i * setm->Nstate_rubrene + i] = sum_val;

            if (i < setm->Nstate_rubrene - 1) {
                Heff[i * setm->Nstate_rubrene + (i + 1)] = setm->Vc_rubrene;
                Heff[(i + 1) * setm->Nstate_rubrene + i] = setm->Vc_rubrene;
            } else {
                Heff[0 * setm->Nstate_rubrene + (setm->Nstate_rubrene - 1)] = setm->Vc_rubrene;
                Heff[(setm->Nstate_rubrene - 1) * setm->Nstate_rubrene + 0] = setm->Vc_rubrene;
            }
        }



        dia_symmat(setm->Nstate_rubrene, Heff, E, C);
 
        memset(expe, 0, setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        for (i = 0; i < setm->Nstate_rubrene; i++) {
            expe[i * setm->Nstate_rubrene + i] = exp(-beta * E[i]);
        }
        tempdm1 = (double *)malloc(setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        tempdm2 = (double *)malloc(setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        transpose(C,tempdm1,setm->Nstate_rubrene);
        dd_matmul(C,expe,tempdm2,setm->Nstate_rubrene,setm->Nstate_rubrene,setm->Nstate_rubrene);
        dd_matmul(tempdm2,tempdm1,rho,setm->Nstate_rubrene,setm->Nstate_rubrene,setm->Nstate_rubrene);
 
        memset(w0, 0, setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        memset(wt, 0, setm->Nstate_rubrene * setm->Nstate_rubrene * sizeof(double));
        
        for (i = 0; i < setm->Nstate_rubrene; i++) {
            if (i == 0) {
                wt[(i + 1) * setm->Nstate_rubrene + i] = 1;
                wt[i * setm->Nstate_rubrene + (i + 1)] = -1;
            } else if (i == setm->Nstate_rubrene - 1) {
                wt[0 * setm->Nstate_rubrene + i] = 1;
                wt[i * setm->Nstate_rubrene + 0] = -1;
            } else {
                wt[(i + 1) * setm->Nstate_rubrene + i] = 1;
                wt[i * setm->Nstate_rubrene + (i + 1)] = -1;
            }
        }
        dd_matmul(rho,wt,w0,setm->Nstate_rubrene,setm->Nstate_rubrene,setm->Nstate_rubrene);

    }

    free(rho);
    free(Heff);
    free(E);
    free(C);
    free(expe);
    free(tempdm1);
    free(tempdm2);

}


