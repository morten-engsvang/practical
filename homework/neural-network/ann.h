#ifndef HAVE_ANN_H
#define HAVE_ANN_H

typedef struct {int n; double (*f)(double x); double (*f_dif)(double x); double (*f_int)(double x); gsl_vector* params; } ann;

ann* ann_alloc(
	int n, /* Antallet af neuroner i hidden layer */
	double f(double x), /* Aktiveringsfunktionen */
	double f_dif(double x), /* Den afledte af aktiveringsfunktionen */
	double f_int(double x) /* Den integrerede af aktiveringsfunktionen */
	);

void ann_free (ann* network);

double ann_response(ann* network, double x);

double ann_response_diff(ann* network, double x);

double ann_response_int(ann* network, double x);

void ann_train(ann* network, gsl_vector* xs, gsl_vector* ys);

#endif
