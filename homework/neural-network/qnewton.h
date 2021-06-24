#ifndef HAVE_QNEWTON_H
#define HAVE_QNEWTON_H

int qnewton(void F(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps);

#endif

