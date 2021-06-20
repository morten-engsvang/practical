#ifndef HAVE_INTEGRATE_H
#define HAVE_INTEGRATE_H

double recursive_integrate(
	double f(double), /* Funktionen der skal integreres */
	double a, /* nederste integrationsgrænse */
	double b, /* øverste integrationsgrænse */
	double delta, /* Tilladte absolutte fejl */
	double eps /* Tilladte relative fejl */
);

double clenshaw_curtis(
	double f(double),
	double a,
	double b,
	double delta,
	double eps
);

#endif
