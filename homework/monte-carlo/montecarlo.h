#ifndef HAVE_MONTECARLO_H
#define HAVE_MONTECARLO_H

void plainmc(
	int dim, /* Antallet af dimensioner */
	double f(int dim,gsl_vector* x), /* Funktionen der skal integreres */
	gsl_vector* a, /* Nederste integrationsgrænser */
	gsl_vector* b, /* Øverste integrationsgrænser */
	int N, /* Antal punkter der skal samples */
	gsl_vector* result_error /* 2D vektor hvor 1. element er resultatet givet i ligning 2, og 2. element er fejlen givet i ligning 3. */
);

#endif
