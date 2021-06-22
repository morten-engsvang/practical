#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include"qr.h"
#include"rkstep12.h"
#include<float.h> //Giver machine epsilon som DBL_EPSILON
double delta = sqrt(DBL_EPSILON);

static double epsilon;

void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps) {
	//I will construct the Jacobian using finite differences, delta_x = sqrt(DBL_EPSILON)
	//And diagonalise it using the routine implemented in a previous homework.
	int n=x->size; //Finder længden af x.
	int maxsteps = 10000;
	int step = 0;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* newx = gsl_vector_alloc(n);
	gsl_vector* backx = gsl_vector_alloc(n);
	gsl_vector* backfx = gsl_vector_alloc(n);
	
	while (step < maxsteps){
		f(x,fx); //Fylder fx med funktionsværdierne.
		//Laver min Jacobian via. finite differences:
		for (int j = 0; j < n; j++){
			double xj = gsl_vector_get(x,j); //Gemmer elementets originalværdi til senere brug.
			gsl_vector_set(x,j,xj+delta);
			f(x,df);
			gsl_vector_sub(df,fx); /* Giver df = f(x+delta)-f(x) for denne værdi af k=j i ligning 7, næste loop over i udvælger så alle ij elementer.*/
			for (int i = 0; i < n; i++){
				gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/delta);
			}
			gsl_vector_set(x,j,xj); //Putter originalværdien tilbage
		}
		GS_decomp(J,R); //Foretager QR faktorisering
		GS_solve(J,R,fx,newx); //Løser JDx=f(x) og placerer løsningen i newx.
		gsl_vector_scale(newx,-1.0); //Løsningen til JDx=-f(x) er den negative af ovenstående løsning.
		//Starter nu backtracking linesearch
		double l = 1;
		while (l > 1./64){ //Backtracker indtil l er 1/64 eller mindre.
			gsl_vector_memcpy(backx,x); //Vi skal bruge x igen og GSL er destruktiv :(
			gsl_vector_add(backx,newx); //laver x+dx vektoren ud fra den tidligere løsning, fuld skridt.
			f(backx,backfx); //Finder funktionsværdien i f(x+dx), dvs. fuld skridt.
			//Tjekker så om den opfylder ligning 8, gsl_blas_dnrm2 giver normen.
			if ( gsl_blas_dnrm2(backfx) < (1-l/2)*gsl_blas_dnrm2(fx)){
				break;
			}
			l *= 0.5; //Vi var gået for langt tilbage, så vi rykker frem.
		}
		gsl_vector_memcpy(backx,x); //Jeg laver kopi af x
		gsl_vector_scale(newx,l); //Jeg skalerer min opdatering til x med den fundne l værdi
		gsl_vector_add(backx,newx); // Lægger min opdatering til x
		gsl_vector_memcpy(x,backx); // Erstatter x med det nye x.
		f(x,fx); //Finder funktionsværdien i det nye punkt.
		//Jeg kan nu tjekke om det har ændret sig nok, tjekker om ændringerne er store nok,
		//eller om vi er tæt nok på til at f(x+dx) = 0
		if ( gsl_blas_dnrm2(newx) < delta || gsl_blas_dnrm2(fx) < eps){
			break;
		}
		step++; //Ej at forglemme
	}
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(df);
	gsl_vector_free(newx);
	gsl_vector_free(backx);
	gsl_vector_free(backfx);
}



}


int main(void){
	
return 0;
}
