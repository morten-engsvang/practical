#ifndef HAVE_SPLINE_H
#define HAVE_SPLINE_H
double linterp(gsl_vector* x, gsl_vector* y, double z, int n);

double linear_integral(double x1, double x2, double y1, double slope);

double linterp_integ(gsl_vector* x, gsl_vector* y, double z, int n);

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline * quad_alloc(int n, double *x, double *y);

void quad_free(qspline *s);

double quad_interp(qspline *s, double z);

double quad_deriv(qspline * s, double z);

double quad_integ(qspline * s, double z);
#endif
