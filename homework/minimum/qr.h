#ifndef HAVE_QR_H
#define HAVE_QR_H

void GS_decomp(
	gsl_matrix* A, /* A er n*m (n>=m) matricen der ønskes dekomponeret til QR, Q erstatter A matricen */
	gsl_matrix* R /* Tom m*m matrice, der erstattes med R matricen.  */
);

void GS_solve(
	gsl_matrix* Q, /* Q matrix fra GS_decomp */
	gsl_matrix* R, /* R matrix fra GS_decomp */
	gsl_vector* b, /* b vektor fra ligningssystemet Ax=b */
	gsl_vector* x /* Tom n vektor, hvor løsningen placeres*/
);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

#endif
