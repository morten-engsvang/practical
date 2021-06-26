#include<gsl/gsl_matrix.h>
#ifndef HAVE_LEAST_SQUARES_H
#define HAVE_LEAST_SQUARES_H
void backsub(gsl_matrix* U, gsl_vector* c);

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

void least_squares(int m, double f(int i,double x), gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c,gsl_matrix* S);

#endif
