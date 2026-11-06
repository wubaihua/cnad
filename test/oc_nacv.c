#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "oc.h"

static int check_case(int n, int forcetype, double alpha, double gamma) {
    struct set_host setm;
    int ndof = 2 * n;
    size_t matrix_size = (size_t)n * n;
    size_t tensor_size = matrix_size * ndof;
    double *R = malloc(ndof * sizeof(*R));
    double *real_dV = malloc(tensor_size * sizeof(*real_dV));
    double *E = malloc(n * sizeof(*E));
    double complex *U = malloc(matrix_size * sizeof(*U));
    double complex *dV = malloc(tensor_size * sizeof(*dV));
    double complex *dv = malloc(tensor_size * sizeof(*dv));
    double complex *nac = malloc(tensor_size * sizeof(*nac));
    double complex *temp = malloc(matrix_size * sizeof(*temp));
    int ok = 1;
    const double pi = acos(-1.0);

    if (!R || !real_dV || !E || !U || !dV || !dv || !nac || !temp) {
        fprintf(stderr, "test allocation failed\n");
        exit(EXIT_FAILURE);
    }

    memset(&setm, 0, sizeof(setm));
    setm.Nmole_oc = n;
    setm.Nstate = n;
    setm.Ndof1 = n;
    setm.Ndof2 = 2;
    setm.forcetype = forcetype;
    setm.alpha_oc = alpha;
    setm.gamma_oc = gamma;
    setm.K_oc[0] = 0.13;
    setm.K_oc[1] = 0.21;

    for (int i = 0; i < ndof; ++i) R[i] = 0.7 * sin(i * 0.37 + 0.1);
    for (int i = 0; i < n; ++i) {
        E[i] = i == 1 ? 1e-8 : 0.05 * i;
        for (int j = 0; j < n; ++j) {
            double phase = 2.0 * pi * i * j / n + 0.13 * j;
            U[i * n + j] = (cos(phase) + I * sin(phase)) / sqrt((double)n);
        }
    }
    dV_OC(R, real_dV, forcetype, &setm);
    for (size_t i = 0; i < tensor_size; ++i) dV[i] = real_dV[i];
    cal_NACV_oc(U, dV, E, dv, nac, &setm);

    for (int k = 0; k < ndof; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double complex sum = 0.0;
                for (int l = 0; l < n; ++l) {
                    sum += conj(U[l * n + i]) * dV[(l * n + j) * ndof + k];
                }
                temp[i * n + j] = sum;
            }
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double complex expected = 0.0;
                size_t index = ((size_t)i * n + j) * ndof + k;
                for (int l = 0; l < n; ++l) expected += temp[i * n + l] * U[l * n + j];
                if (cabs(dv[index] - expected) > 2e-12 * (1.0 + cabs(expected)) ||
                    (i == j && cabs(nac[index]) != 0.0) ||
                    (i != j && cabs((nac[index] - expected / (E[j] - E[i]))
                                    * (E[j] - E[i])) > 2e-12 * (1.0 + cabs(expected)))) {
                    fprintf(stderr, "n=%d force=%d alpha=%g gamma=%g i=%d j=%d k=%d mismatch\n",
                            n, forcetype, alpha, gamma, i, j, k);
                    ok = 0;
                    goto done;
                }
            }
        }
    }

done:
    free(R);
    free(real_dV);
    free(E);
    free(U);
    free(dV);
    free(dv);
    free(nac);
    free(temp);
    return ok;
}

int main(void) {
    int ok = 1;
    const int sizes[] = {3, 7, 51};
    for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s) {
        for (int forcetype = 0; forcetype <= 1; ++forcetype) {
            ok &= check_case(sizes[s], forcetype, 0.8, 0.3);
            ok &= check_case(sizes[s], forcetype, 0.0, 0.0);
        }
    }
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
