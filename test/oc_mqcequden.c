#include <math.h>
#include <stdio.h>
#include <string.h>

#include "oc.h"

static int close_enough(double actual, double expected) {
    return isfinite(actual) && fabs(actual - expected) < 1e-11 * (1.0 + fabs(expected));
}

static int check_matrix(const char *name, const double *actual, const double *expected) {
    for (int i = 0; i < 9; ++i) {
        if (!close_enough(actual[i], expected[i])) {
            fprintf(stderr, "%s[%d]: got %.17g, expected %.17g\n",
                    name, i, actual[i], expected[i]);
            return 0;
        }
    }
    return 1;
}

int main(void) {
    struct set_host setm;
    double R[6] = {0.4, 0.3, -0.2, -0.4, 0.7, 0.8};
    double P[6] = {0.8, -0.5, 0.4, 0.2, -0.6, 0.9};
    double rho[9], expected[9], electronic[9];
    double beta = 0.7;
    int ok = 1;

    memset(&setm, 0, sizeof(setm));
    setm.Nmole_oc = 3;
    setm.mass_oc[0] = 2.0;
    setm.mass_oc[1] = 4.0;
    setm.K_oc[0] = 1.5;
    setm.K_oc[1] = 0.8;
    setm.alpha_oc = 0.6;
    setm.if_mqcequden = 2;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            expected[i * 3 + j] = i == j ? exp(-beta * setm.alpha_oc * R[2 * i]) : 0.0;
        }
    }
    mqcequden_OC(rho, beta, R, P, &setm);
    ok &= check_matrix("diagonal electronic", rho, expected);

    memcpy(electronic, rho, sizeof(rho));
    double nuclear_energy = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int mode = 0; mode < 2; ++mode) {
            int coordinate = 2 * i + mode;
            nuclear_energy += 0.5 * (setm.K_oc[mode] * R[coordinate] * R[coordinate]
                                     + P[coordinate] * P[coordinate] / setm.mass_oc[mode]);
        }
    }
    for (int i = 0; i < 9; ++i) {
        expected[i] = electronic[i] * exp(-beta * nuclear_energy);
    }
    setm.if_mqcequden = 1;
    mqcequden_OC(rho, beta, R, P, &setm);
    ok &= check_matrix("full phase-space", rho, expected);

    setm.if_mqcequden = 2;
    setm.alpha_oc = 0.0;
    setm.beta_oc = 0.4;
    double ground_weight = exp(-2.0 * beta * setm.beta_oc);
    double degenerate_weight = exp(beta * setm.beta_oc);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            expected[i * 3 + j] = (ground_weight - degenerate_weight) / 3.0
                                + (i == j ? degenerate_weight : 0.0);
        }
    }
    mqcequden_OC(rho, beta, R, P, &setm);
    ok &= check_matrix("coupled ring", rho, expected);

    setm.alpha_oc = 0.6;
    setm.gamma_oc = 0.3;
    mqcequden_OC(rho, beta, R, P, &setm);
    {
        double H[9];
        V_OC(R, H, 1, &setm);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double rh = 0.0, hr = 0.0;
                for (int k = 0; k < 3; ++k) {
                    rh += rho[i * 3 + k] * H[k * 3 + j];
                    hr += H[i * 3 + k] * rho[k * 3 + j];
                }
                if (!close_enough(rho[i * 3 + j], rho[j * 3 + i]) ||
                    !close_enough(rh, hr)) {
                    fprintf(stderr, "inhomogeneous ring matrix check failed at %d,%d\n", i, j);
                    ok = 0;
                }
            }
        }
    }

    for (int i = 0; i < 9; ++i) {
        expected[i] = i % 4 == 0 ? 1.0 : 0.0;
    }
    mqcequden_OC(rho, 0.0, R, P, &setm);
    ok &= check_matrix("zero beta", rho, expected);
    return ok ? 0 : 1;
}
