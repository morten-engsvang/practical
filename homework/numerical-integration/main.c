#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include"integrate.h"

double testf1(double x){
	return sqrt(x);
}

double testf2(double x){
	double U = 4*sqrt(1-pow(x,2));
	return U;
}

double testf3(double x){
	double U = 1.0/sqrt(x);
	return U;
}

double testf4(double x){
	double U = log(x)/sqrt(x);
	return U;
}



int main(void){
	double delta = 0.01;
	double eps = 0.01;
	int evals = 0;
	double a = 0;
	double b = 1;
	/* Jeg laver min tests til del A */
	double function1(double x){
		evals++;
		return testf1(x);
	}
	printf("Test 1: f(x) = sqrt(x) from %g to %g without Clenshaw-Curtis\n",a,b);
	double test1 = recursive_integrate(function1,a,b,delta,eps);
	printf("The result of test 1: %g\n",test1);
	printf("The number of evaluations in test 1: %i\n",evals);
	printf("Test 2: f(x) = 4*sqrt(1-x^2) from %g to %g without Clenshaw-Curtis\n",a,b);
	evals = 0;
	double function2(double x){
		evals++;
		return testf2(x);
	}
	double test2 = recursive_integrate(function2,a,b,delta,eps);
	printf("The result of test 2: %g\n",test2);
	printf("The number of evaluations in test 2: %i\n",evals);
	
	/* Jeg laver mine tests til del b */
	//første integral uden clenshaw
	
	
	printf("Test 3: f(x) = 1/sqrt(x) from %g to %g without Clenshaw-Curtis\n",a,b);
	evals = 0;
	double function3(double x){
		evals++;
		return testf3(x);
	}
	double test3 = recursive_integrate(function3,a,b,delta,eps);
	printf("The result of test 3: %g\n",test3);
	printf("The number of evaluations in test 3: %i\n",evals);
	//første integral med clenshaw
	printf("Test 4: f(x) = 1/sqrt(x) from %g to %g with Clenshaw-Curtis\n",a,b);
	evals = 0;
	
	
	/*
	double test4 = clenshaw_curtis(function3,a,b,delta,eps);
	printf("The result of test 4: %g\n",test4);
	printf("The number of evaluations in test 4: %i\n",evals);
	*/
	
	//andet integral uden clenshaw
	
	/*
	printf("Test 5: f(x) = ln(x)/sqrt(x) from %g to %g without Clenshaw-Curtis\n",a,b);
	evals = 0;
	double function4(double x){
		evals++;
		return testf4(x);
	}
	double test5 = recursive_integrate(function4,a,b,delta,eps);
	printf("The result of test 5: %g\n",test5);
	printf("The number of evaluations in test 3: %i\n",evals);
	*/
	
	//andet integral med clenshaw
	
	/*
	printf("Test 6: f(x) = 1/sqrt(x) from %g to %g with Clenshaw-Curtis\n",a,b);
	evals = 0;
	double test6 = clenshaw_curtis(function4,a,b,delta,eps);
	printf("The result of test 6: %g\n",test6);
	printf("The number of evaluations in test 4: %i\n",evals);
	*/
	
return 0;
}
