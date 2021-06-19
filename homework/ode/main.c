#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include"rkstep12.h"
void test_function(double t, gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}

void SIR(double t, gsl_vector* y, gsl_vector* dydt){
	double TC = 3;
	double TR = 20;
	double N;
	double S = gsl_vector_get(y,0); 
	double I = gsl_vector_get(y,1); 
	double R = gsl_vector_get(y,2);
	N = S + I + R;
	gsl_vector_set(dydt,0,-I*S/(N*TC));
	gsl_vector_set(dydt,1,I*S/(N*TC)-I/TR);
	gsl_vector_set(dydt,2,I/TR);
}


int main(void){
	/* Bruger den givne testfunktion: */
	double a = 0;
	double b = 10;
	double h = 0.02;
	int n = 2;
	double acc = 0.001;
	double eps = 0.001;
	FILE* testout = fopen("test.txt","w");
	gsl_vector* ya = gsl_vector_alloc(n);
	gsl_vector* yb = gsl_vector_alloc(n);
	gsl_vector_set(ya,0,1); gsl_vector_set(ya,1,0);
	//printf("Klar til at køre driver\n");
	driver(test_function, a, ya, b, yb, h, acc, eps, n, testout);
	//printf("The final value of f(b): %g\n", gsl_vector_get(yb,0));
	//printf("The final value of f'(b): %g\n", gsl_vector_get(yb,1));
	fclose(testout);
	
	gsl_vector_free(ya);
	gsl_vector_free(yb);
	
	/* Jeg kan nu teste SIR modellen */
	a = 0;
	b = 100;
	h = 0.02;
	n = 3;
	acc = 0.001;
	eps = 0.001;
	gsl_vector* ystart = gsl_vector_alloc(n);
	gsl_vector* yslut = gsl_vector_alloc(n);
	gsl_vector_set(ystart,0,1000000);
	gsl_vector_set(ystart,1,1000);
	gsl_vector_set(ystart,2,0);
	FILE* sirout = fopen("sir.txt","w");
	driver(SIR,a,ystart,b,yslut,h,acc,eps,n,sirout);
	fclose(sirout);
	gsl_vector_free(ystart);
	gsl_vector_free(yslut);
return 0;
}

