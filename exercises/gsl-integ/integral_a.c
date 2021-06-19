#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

//Jeg definerer funktionen der ønskes integreret

double function (double x, void* params) {
	double f = log(x)/sqrt(x);
	return f;
}


int main(void) {
	double z;
	gsl_function F;
	F.function=&function;
	F.params=(void*)&z;
	int limit = 999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0,b=1,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	printf("The approximated value of integral a is: %g\n",result);
return 0;
}
