#ifndef HAVE_JACOBI_H
#define HAVE_JACOBI_H

void jacobi_diag(
	gsl_matrix* A, /* Matrix der ønskes diagonaliseret */
	gsl_matrix* V /* Enhedsmatrice */
	/* Output: A er diagonaliseret in-place, V er nu en matrix af egenvektorer. */
);

#endif
