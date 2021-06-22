#ifndef HAVE_RKSTEP_H
#define HAVE_RKSTEP_H
void rkstep12(
	void f(double t, gsl_vector* y, gsl_vector* dydt), /* The f from dy/dt = f(t,y) */
	double t, /* Current value of t */
	gsl_vector* yt, /*  The current value y(t) of the sought function */
	double h, /*  The step to be taken */
	gsl_vector* yh, /*  output: y(t+h)  */
	gsl_vector* err, /*  output: error estimate  */
	int n /* Antallet af koblede ligninger */
);

void driver(
	void f(double t, gsl_vector* y, gsl_vector* dydt), /* right-hand-side of dy/dt=f(t,y) */
	double a, /* the start-point a */
	gsl_vector* ya, /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector* yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,                    /* relative accuracy goal */
	int n, /* Antallet af koblede ligninger */
	FILE * outstream
);

#endif
