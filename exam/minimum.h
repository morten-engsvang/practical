#include<gsl/gsl_matrix.h>
#ifndef HAVE_MINIMUM_H
#define HAVE_MINIMUM_H
int qnewton(void F(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps);

#endif
