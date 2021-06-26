#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include"montecarlo.h"

double testf1(int dim, gsl_vector* x){
	double result = 1.0;
	for (int i = 0; i < dim; i++){
		double U = 4*sqrt(1-pow(gsl_vector_get(x,i),2));
		result *= U;
	}
	
	//double result = 4*sqrt(1-pow(gsl_vector_get(x,dim),2));
	
	return result;
	//Expected result when integrated from a to b in 1D is pi.
}

double integral(int dim, gsl_vector* x){
	double temp = 1-cos(gsl_vector_get(x,0))*cos(gsl_vector_get(x,1))*cos(gsl_vector_get(x,2));
	double result = pow(temp,-1.0)/pow(M_PI,3);
	return result;
}


int main(void){
	//Tester først om den kan beregne i 1D
	printf("This is part A:\n");
	printf("Test 1 is: U(x) = 4*sqrt(1-x^2) from 0 to 1, to see if it works at all.\n");
	int dim = 1;
	gsl_vector* a = gsl_vector_alloc(1);
	gsl_vector* b = gsl_vector_alloc(1);
	gsl_vector* result_error = gsl_vector_alloc(2);
	gsl_vector_set(a,0,0);
	gsl_vector_set(b,0,1);
	int N = 100000;
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
	
	//Beregner det integral han gerne vil have.
	printf("Test 3 is the integral given in the homework.\n");
	gsl_vector* a3 = gsl_vector_alloc(dim);
	gsl_vector* b3 = gsl_vector_alloc(dim);
	dim = 3;
	long int N1 = 10000000;
	for (int i = 0; i < dim; i++){
		gsl_vector_set(a3,i,0);
		gsl_vector_set(b3,i,M_PI);
	}
	plainmc(dim,integral,a3,b3,N1,result_error);
	printf("The result of test 3 is: %g\n",gsl_vector_get(result_error,0));
	printf("The estimated error of test 3 is: %g\n",gsl_vector_get(result_error,1));
	printf("I used %li points.\n",N1);
	
	printf("-----------------------\n");
	printf("This is part B\n");
	printf("I estimate error using two different sequences: err = fabs(result1-result2)/(N*volume)\n");
	printf("where result1 and result2 is from sampling N/2 points each for a total of N points.\n");
	printf("I can now test the scaling of the error, I have plotted the error in plot.svg and plot_actual.svg\n");
	printf("Where I have plotted the error estimates and the actual errors respectively\n");
	printf("I test it for the integral given in the homework\n");
	
	double actual = 1.3932039296856768591842462603255;
	printf("Which means that the actual error is relative to: %g",actual);
	dim = 3;
	FILE * data = fopen("plot.txt","w");
	for (int N2 = 10; N2 < 500; N2 += 5){
		plainmc(dim,integral,a3,b3,N2,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral,a3,b3,N2,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data, "%i %10g %10g %10g %10g\n", N2, p_est_error, p_act_error, h_est_error, h_act_error);
	}
	fclose(data);
	printf("It can seen that the plain Monte-Carlo integration error falls relatively quickly\nafter which it begins oscillating around 0\n");
	printf("The error using the error using the quasi-random Halton sequence starts much lower and\nappears to be more stable. However it has weird dips in the actual error.\n");
	gsl_vector_free(a);
	gsl_vector_free(b);
	gsl_vector_free(result_error);
	gsl_vector_free(a2);
	gsl_vector_free(b2);
	gsl_vector_free(a3);
	gsl_vector_free(b3);
return 0;
}
