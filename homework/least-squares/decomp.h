#ifndef HAVE_DECOMP_H
#define HAVE_DECOMP_H
void backsub(gsl_matrix* U, gsl_vector* c);

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

#endif
