#ifndef GMATH_H
#define GMATH_H

#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "constant.h"  // Assuming constants are defined in a separate header
#include <lapacke.h>

// Generate two random numbers x1, x2 with Gaussian distribution N(miu, sigma)
void box_muller(double *x1, double *x2, double sigma, double miu);

// Kronecker delta function
int Kronecker_delta(int i, int j);

// Heaviside step function
int heaviside(double x);
// Generate permutation matrix U
void gen_permutation(int n, double **U, int offset);

// Sort array a of length n
void math_sort(int n, double *a);
// Generate a random probability distribution p of length n
void random_prob(int n, double *p);

// Diagonalize a symmetric matrix A (real, n x n)
void dia_symmat(int n, double *A, double *E, double *C);

// Diagonalize a Hermitian matrix A (complex, n x n)
void dia_hermitemat(int n, double complex *A, double *E, double complex *C);

// Trace of a square matrix A (real, n x n)
double trace(int n, double *A) ;

// Trace of a square matrix A (complex, n x n)
double complex trace_comp(int n, double complex *A) ;

void init_seed(int my_prl);


void dd_matmul(double *A, double *B, double *C, int nA, int nB, int nC);




void cc_matmul(double complex *A, double complex *B, double complex *C, int nA, int nB, int nC);



void dc_matmul(double *A, double complex *B, double complex *C, int nA, int nB, int nC);
void cd_matmul(double complex *A, double *B, double complex *C, int nA, int nB, int nC);

void transpose(double *A, double *AT, int n) ;
void transpose_conjugate(double complex *A, double complex *AT, int n);
void diagger(double complex *A, double complex *Ad, int n) ;

int any_isnan(double complex *array, int size) ;

int maxloc(double *array, int size);

// double complex cc_multiply(double omplex a, double complex b);

// double complex dc_multiply(double real, double complex c);

#endif // GMATH_H
