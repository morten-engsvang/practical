#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include"montecarlo.h"

double testf1(int dim, gsl_vector* x){
	/*
	double result = 1.0;
	for (int i = 0; i < dim; i++){
		double U = 4*sqrt(1-pow(gsl_vector_get(x,dim-1),2));
		result *= U;
	}
	*/
	double result = 4*sqrt(1-pow(gsl_vector_get(x,dim),2));
	
	return result;
	//Expected result when integrated from a to b in 1D is pi.
}


int main(void){
	//Tester først om den kan beregne i 1D
	printf("Test 1 is: U(x) = 4*sqrt(1-x^2) from 0 to 1, to see if it works at all.\n");
	int dim = 1;
	gsl_vector* a = gsl_vector_alloc(1);
	gsl_vector* b = gsl_vector_alloc(1);
	gsl_vector* result_error = gsl_vector_alloc(2);
	gsl_vector_set(a,0,0);
	gsl_vector_set(b,0,1);
	int N = 1000;
	plainmc(dim,testf1,a,b,N,result_error);
	printf("The result of test 1 is: %g\n",gsl_vector_get(result_error,0));
	printf("The estimated error of test 1 is: %g\n",gsl_vector_get(result_error,1));
	printf("The expected result is: %g and the actual error is then: %g\n",M_PI,fabs(gsl_vector_get(result_error,0)-M_PI));
	
	
	//Tester så lige om den virker i flere dimensioner:
	dim = 3;
	gsl_vector* a2 = gsl_vector_alloc(dim);
	gsl_vector* b2 = gsl_vector_alloc(dim);
	for (int i = 0; i < dim; i++){
		gsl_vector_set(a2,i,0);
		gsl_vector_set(b2,i,1);
	}
	printf("Test 2 is the same U(x) but now f(x,y,z) = U(x)*U(y)*U(z)\n");
	plainmc(dim,testf1,a2,b2,N,result_error);
	printf("The result of test 2 is: %g\n",gsl_vector_get(result_error,0));
	printf("The estimated error of test 2 is: %g\n",gsl_vector_get(result_error,1));
	printf("The expected result is: %g and the actual error is then: %g\n",pow(M_PI,dim),fabs(gsl_vector_get(result_error,0)-pow(M_PI,dim)));
	
	
	gsl_vector_free(a);
	gsl_vector_free(b);
	gsl_vector_free(result_error);
	gsl_vector_free(a2);
	gsl_vector_free(b2);
return 0;
}
