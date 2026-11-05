#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constant.h"
#include "gmath.h"
#include "oc.h"

/*
 * Wang and Beljonne, J. Phys. Chem. Lett. 2013, 4, 1888-1894:
 *
 * H_e = sum_k alpha x_(k,1) |k><k|
 *     + sum_k [beta - gamma (x_(k+1,2) - x_(k,2))]
 *       (|k><k+1| + |k+1><k|)
 *
 * H_n = sum_(k,s) [m_s v_(k,s)^2 + K_s x_(k,s)^2] / 2.
 *
 * Coordinates are stored in bohr and all parameters in setm are converted to
 * atomic units by readinp_OC.
 */

static int oc_coordinate(int site, int mode) {
    return 2 * site + mode;
}

void parameter_OC(double *mass, struct set_host *setm) {
    int site;

    for (site = 0; site < setm->Nmole_oc; ++site) {
        mass[oc_coordinate(site, 0)] = setm->mass_oc[0];
        mass[oc_coordinate(site, 1)] = setm->mass_oc[1];
    }
}

void sample_OC(double *P, double *R, double beta, struct set_host *setm) {
    int site, mode;
    double x2;

    for (site = 0; site < setm->Nmole_oc; ++site) {
        for (mode = 0; mode < 2; ++mode) {
            int coordinate = oc_coordinate(site, mode);
            double omega = sqrt(setm->K_oc[mode] / setm->mass_oc[mode]);

            if (setm->if_classical) {
                box_muller(&P[coordinate], &x2, sqrt(setm->mass_oc[mode] / beta), 0.0);
                box_muller(&R[coordinate], &x2, sqrt(1.0 / (beta * setm->K_oc[mode])), 0.0);
            } else if (beta > 99999.0) {
                box_muller(&P[coordinate], &x2,
                           sqrt(0.5 * hbar * setm->mass_oc[mode] * omega), 0.0);
                box_muller(&R[coordinate], &x2,
                           sqrt(0.5 * hbar / (setm->mass_oc[mode] * omega)), 0.0);
            } else {
                double thermal_factor = tanh(0.5 * beta * hbar * omega);

                box_muller(&P[coordinate], &x2,
                           sqrt(0.5 * hbar * setm->mass_oc[mode] * omega / thermal_factor), 0.0);
                box_muller(&R[coordinate], &x2,
                           sqrt(0.5 * hbar / (thermal_factor * setm->mass_oc[mode] * omega)), 0.0);
            }
        }
    }
}

void V_OC(double *R, double *H, int forcetype, struct set_host *setm) {
    int i, j;
    int nsite = setm->Nmole_oc;

    memset(H, 0, nsite * nsite * sizeof(double));

    for (i = 0; i < nsite; ++i) {
        int next = (i + 1) % nsite;
        double coupling = setm->beta_oc
                        - setm->gamma_oc
                        * (R[oc_coordinate(next, 1)] - R[oc_coordinate(i, 1)]);

        H[i * nsite + i] = setm->alpha_oc * R[oc_coordinate(i, 0)];
        H[i * nsite + next] += coupling;
        H[next * nsite + i] += coupling;
    }

    if (forcetype == 0) {
        double nuclear_potential = 0.0;

        for (i = 0; i < nsite; ++i) {
            for (j = 0; j < 2; ++j) {
                double coordinate = R[oc_coordinate(i, j)];
                nuclear_potential += 0.5 * setm->K_oc[j] * coordinate * coordinate;
            }
        }

        for (i = 0; i < nsite; ++i) {
            H[i * nsite + i] += nuclear_potential;
        }
    }
}

void dV_OC(double *R, double *dH, int forcetype, struct set_host *setm) {
    int i, j;
    int nsite = setm->Nmole_oc;
    int ndof = 2 * nsite;

    memset(dH, 0, nsite * nsite * ndof * sizeof(double));

    for (i = 0; i < nsite; ++i) {
        int next = (i + 1) % nsite;
        int local_coordinate = oc_coordinate(i, 0);
        int nonlocal_coordinate = oc_coordinate(i, 1);
        int next_nonlocal_coordinate = oc_coordinate(next, 1);

        dH[(i * nsite + i) * ndof + local_coordinate] = setm->alpha_oc;

        dH[(i * nsite + next) * ndof + nonlocal_coordinate] += setm->gamma_oc;
        dH[(next * nsite + i) * ndof + nonlocal_coordinate] += setm->gamma_oc;
        dH[(i * nsite + next) * ndof + next_nonlocal_coordinate] -= setm->gamma_oc;
        dH[(next * nsite + i) * ndof + next_nonlocal_coordinate] -= setm->gamma_oc;
    }

    if (forcetype == 0) {
        for (i = 0; i < nsite; ++i) {
            for (j = 0; j < nsite; ++j) {
                int mode;

                for (mode = 0; mode < 2; ++mode) {
                    int coordinate = oc_coordinate(i, mode);
                    dH[(j * nsite + j) * ndof + coordinate] += setm->K_oc[mode] * R[coordinate];
                }
            }
        }
    }
}

void nucforce_OC(double *R, double *nf, struct set_host *setm) {
    int site, mode;

    for (site = 0; site < setm->Nmole_oc; ++site) {
        for (mode = 0; mode < 2; ++mode) {
            int coordinate = oc_coordinate(site, mode);
            nf[coordinate] = setm->K_oc[mode] * R[coordinate];
        }
    }
}

static void oc_current_operator(double *R, double *current, struct set_host *setm) {
    int site;
    int nsite = setm->Nmole_oc;

    memset(current, 0, nsite * nsite * sizeof(double));
    for (site = 0; site < nsite; ++site) {
        int next = (site + 1) % nsite;
        double coupling = setm->beta_oc
                        - setm->gamma_oc
                        * (R[oc_coordinate(next, 1)] - R[oc_coordinate(site, 1)]);

        current[next * nsite + site] += coupling;
        current[site * nsite + next] -= coupling;
    }
}

void cfweight_OC(double *w0, double *wt, double beta, double *R, double *P,
                 double *Rt, struct set_host *setm) {
    int i;
    int nsite = setm->Nmole_oc;
    double *rho = (double *)malloc(nsite * nsite * sizeof(double));
    double *Heff = (double *)malloc(nsite * nsite * sizeof(double));
    double *E = (double *)malloc(nsite * sizeof(double));
    double *C = (double *)malloc(nsite * nsite * sizeof(double));
    double *expe = (double *)malloc(nsite * nsite * sizeof(double));
    double *tempdm1 = (double *)malloc(nsite * nsite * sizeof(double));
    double *tempdm2 = (double *)malloc(nsite * nsite * sizeof(double));

    (void)P;

    if (setm->if_classical == 0 || setm->if_classical == 2) {
        double reorganization_energy = 0.5 * setm->alpha_oc * setm->alpha_oc
                                     / setm->K_oc[0];
        double effective_coupling = setm->beta_oc
                                  * exp(-reorganization_energy * beta / 3.0);

        memset(Heff, 0, nsite * nsite * sizeof(double));
        memset(wt, 0, nsite * nsite * sizeof(double));
        for (i = 0; i < nsite; ++i) {
            int next = (i + 1) % nsite;

            Heff[i * nsite + next] += effective_coupling;
            Heff[next * nsite + i] += effective_coupling;
            wt[next * nsite + i] += effective_coupling;
            wt[i * nsite + next] -= effective_coupling;
        }
    } else {
        V_OC(R, Heff, 1, setm);
        oc_current_operator(Rt, wt, setm);
    }

    dia_symmat(nsite, Heff, E, C);

    memset(expe, 0, nsite * nsite * sizeof(double));
    for (i = 0; i < nsite; ++i) {
        expe[i * nsite + i] = exp(-beta * E[i]);
    }

    transpose(C, tempdm1, nsite);
    dd_matmul(C, expe, tempdm2, nsite, nsite, nsite);
    dd_matmul(tempdm2, tempdm1, rho, nsite, nsite, nsite);

    if (setm->if_classical == 0 || setm->if_classical == 2) {
        double trace = 0.0;

        for (i = 0; i < nsite; ++i) {
            trace += rho[i * nsite + i];
        }
        for (i = 0; i < nsite * nsite; ++i) {
            rho[i] /= trace;
        }

        memcpy(Heff, wt, nsite * nsite * sizeof(double));
    } else {
        oc_current_operator(R, Heff, setm);
    }

    dd_matmul(rho, Heff, w0, nsite, nsite, nsite);

    free(rho);
    free(Heff);
    free(E);
    free(C);
    free(expe);
    free(tempdm1);
    free(tempdm2);
}