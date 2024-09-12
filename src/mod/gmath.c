

#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "constant.h"  // Assuming constants are defined in a separate header
#include <lapacke.h>
#include "gmath.h"
#include <string.h>
#include <ctype.h>
#include <stdio.h>

// Generate two random numbers x1, x2 with Gaussian distribution N(miu, sigma)
void box_muller(double *x1, double *x2, double sigma, double miu) {
    double u1, u2;
    u1 = rand() / (double)RAND_MAX;
    u2 = rand() / (double)RAND_MAX;
    *x1 = sqrt(-2 * log(u1)) * cos(2 * pi * u2) * sigma + miu;
    *x2 = sqrt(-2 * log(u1)) * sin(2 * pi * u2) * sigma + miu;
}

// Kronecker delta function
int Kronecker_delta(int i, int j) {
    return (i == j) ? 1 : 0;
}

// Heaviside step function
int heaviside(double x) {
    return (x <= 0) ? 0 : 1;
}

// // Generate permutation matrix U
// void gen_permutation(int n, double **U, int offset) {
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             U[i][j] = 0;
//         }
//     }
//     for (int i = 0; i < n; i++) {
//         int j = (i + offset) % n;
//         if (j < 0) j += n;
//         U[i][j] = 1.0;
//     }
// }

// Sort array a of length n
void math_sort(int n, double *a) {
    double *b = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        int max_idx = 0;
        for (int j = 1; j < n; j++) {
            if (a[j] > a[max_idx]) max_idx = j;
        }
        b[n - i - 1] = a[max_idx];
        a[max_idx] = -INFINITY;  // Mark this as used
    }
    for (int i = 0; i < n; i++) {
        a[i] = b[i];
    }
    free(b);
}

// Generate a random probability distribution p of length n
void random_prob(int n, double *p) {
    double *p0 = (double *)malloc((n - 1) * sizeof(double));
    for (int i = 0; i < n - 1; i++) {
        p0[i] = rand() / (double)RAND_MAX;
    }
    math_sort(n - 1, p0);
    for (int j = 0; j < n; j++) {
        if (j == 0) {
            p[j] = p0[0];
        } else if (j == n - 1) {
            p[j] = 1 - p0[n - 2];
        } else {
            p[j] = p0[j] - p0[j - 1];
        }
    }
    free(p0);
}

// Diagonalize a symmetric matrix A (real, n x n)
void dia_symmat(int n, double *A, double *E, double *C) {
    // double *work = (double *)malloc(3 * n * sizeof(double));
    int n2 = n*n;
    // for (int i = 0; i < n2; i++){
    //     C[i]=A[i];
    // }
    memcpy(C, A, n2*sizeof(double));
    LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', n, C, n, E);
    // free(work);
}

// Diagonalize a Hermitian matrix A (complex, n x n)
void dia_hermitemat(int n, double complex *A, double *E, double complex *C) {
    // double *rwork = (double *)malloc(3 * n * sizeof(double));
    // double complex *work = (double complex *)malloc(3 * n * sizeof(double complex));
    int n2 = n*n;
    memcpy(C, A, n2*sizeof(double complex));
    LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'L', n, C, n, E);
    // free(rwork);
    // free(work);
}

// Trace of a square matrix A (real, n x n)
double trace(int n, double *A) {
    double tA = 0;
    for (int i = 0; i < n; i++) {
        tA += A[i*n+i];
    }
    return tA;
}

// Trace of a square matrix A (complex, n x n)
double complex trace_comp(int n, double complex *A) {
    double complex tA = 0 + 0 * I;
    for (int i = 0; i < n; i++) {
        tA += A[i*n+i];
    }
    return tA;
}

// Initialize random seed
// void init_seed(int my_prl) {
//     int n, ival[8], v[3], i;
//     int *seed;
//     time_t t = time(NULL);
//     struct tm *tm_info = localtime(&t);

//     ival[0] = tm_info->tm_sec;
//     ival[1] = tm_info->tm_min;
//     ival[2] = tm_info->tm_hour;
//     ival[3] = tm_info->tm_mday;
//     ival[4] = tm_info->tm_mon + 1;
//     ival[5] = tm_info->tm_year + 1900;
//     ival[6] = tm_info->tm_wday;
//     ival[7] = tm_info->tm_yday;

//     v[0] = 101 * ival[7] + 256 * ival[6] + ival[7] % 103;
//     v[1] = ival[5] + 64 * ival[4] + (1993 % (2 + ival[7]) + 3) * 97 * ival[7];
//     v[2] = ival[2] + 4 * ival[1] + 16 * ival[0] + (997 % (1 + ival[7]) + 5) * 101 * ival[7];

//     seed = (int *)malloc(n * sizeof(int));
//     srand(time(NULL));

//     for (i = 0; i < n; i++) {
//         seed[i] = rand() + v[i % 3] + ival[7];
//     }

//     if (my_prl) {
//         for (i = 0; i < n; i++) {
//             seed[i] += 113 * my_prl;
//         }
//     }

//     free(seed);
// }


void dd_matmul(double *A, double *B, double *C, int nA, int nB, int nC){
    
    double sum;
    for(int i=0; i < nA; i++){
        for(int k=0; k<nC; k++){
            sum=0;
            for(int j=0; j<nB; j++){
                sum+=A[i*nB+j]*B[j*nC+k];
            }
            C[i*nC+k]=sum;
        }
    }
}




void cc_matmul(double complex *A, double complex *B, double complex *C, int nA, int nB, int nC){
    // if (A == NULL || B == NULL || C == NULL) {
    //     fprintf(stderr, "Error: Null pointer passed to cc_matmul\n");
    //     exit(EXIT_FAILURE);
    // }

    double complex sum;
    for(int i=0; i<nA; i++){
        for(int k=0; k<nC; k++){
            sum=0;
            for(int j=0; j<nB; j++){
                sum+=A[i*nB+j]*B[j*nC+k];
            }
            C[i*nC+k]=sum;

        }
    }
}



void dc_matmul(double *A, double complex *B, double complex *C, int nA, int nB, int nC){
    double complex sum;
    for(int i=0; i<nA; i++){
        for(int k=0; k<nC; k++){
            sum=0;
            for(int j=0; j<nB; j++){
                sum+=A[i*nB+j]*B[j*nC+k];
            }
            C[i*nC+k]=sum;
        }
    }
}

void cd_matmul(double complex *A, double *B, double complex *C, int nA, int nB, int nC){
    double complex sum;
    for(int i=0; i<nA; i++){
        for(int k=0; k<nC; k++){
            sum=0;
            for(int j=0; j<nB; j++){
                sum+=A[i*nB+j]*B[j*nC+k];
            }
            C[i*nC+k]=sum;
        }
    }
}


void transpose(double *A, double *AT, int n) {
   
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            AT[j * n + i] = A[i * n + j];
        }
    }

}



void diagger(double complex *A, double complex *Ad, int n) {
   
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Ad[j * n + i] = conj(A[i * n + j]);
        }
    }

}



int any_isnan(double complex *array, int size) {
    for (int i = 0; i < size; i++) {
        if (isnan(creal(array[i])) || isnan(cimag(array[i]))) {
            return 1;
        }
    }
    return 0;
}



char* trim(char *str) {
    char *end;
    while (isspace((unsigned char)*str)) str++;
    if (*str == 0) return str;
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return str;
}

char* adjustl(char *str) {
    char *start = str;
    while (*start == ' ') start++;
    return start;
}

