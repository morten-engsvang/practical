#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include"integrate.h"
#include <gsl/gsl_integration.h>

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



int main(void){
	double delta = 0.005;
	double eps = 0.005;
	int evals = 0;
	double a = 0.0;
	double b = 1.0;
	printf("I have set the maximum allowed amount of recursions at 1000\n");
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
	printf("The error is: %g\n",fabs(test2-M_PI));
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
	double test4 = clenshaw_curtis(function3,a,b,delta,eps);
	printf("The result of test 4: %g\n",test4);
	printf("The number of evaluations in test 4: %i\n",evals);
	//Jeg kan beregne pi igen men nu med Clenshaw-Curtis
	printf("Test 5: repeat of test 2 but with Clenshaw-Curtis\n");
	evals = 0;
	double test5 = clenshaw_curtis(function2,a,b,delta,eps);
	printf("The result of test 5: %g\n",test5);
	printf("The error is: %g\n",fabs(test5-M_PI));
	printf("The number of evaluations in test 5: %i\n",evals);
	//Jeg kan så bruge GSL subroutinen:
	printf("I can calculate the same using the QACS subsroutine from the GSL library\n");
	int limits = 1e6; double result, error; evals = 0;
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(limits);
	gsl_function F;
	double function4(double x, void * params){
		evals++;
		return testf2(x);
	}
	F.function = &function4;
	gsl_integration_qags(&F,a,b,delta,eps,limits,workspace,&result,&error);
	printf("Result, error and number of function evaluations are: \n%.15g\n%.15g\n%i\n",result,fabs(result-M_PI),evals);
	gsl_integration_workspace_free (workspace);
	printf("-------------------------\n");
	printf("A remarkable reduction in evaluations is observed for the integral 1/sqrt(x) from %g to %g\n",a,b);
	printf("Using Clenshaw-Curtis to the pi integral gives a small reduction in error,\nbut also increases the amount of evaluations significantly.\n");
	printf("The GSL routine beats me by a long shot both in error and function evaluations");
return 0;
}
