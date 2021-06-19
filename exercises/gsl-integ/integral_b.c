#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

//Jeg definerer errorfunctionen

double function (double x, void* params) {
	double f = 2/sqrt(M_PI)*exp(-x*x);
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
	double a=0,acc=1e-6,eps=1e-6,result,error;
	double xmin=-5,xmax=5;
	for(double x=xmin;x<=xmax;x+=1.0/8) { 
		gsl_integration_qags(&F,a,x,acc,eps,limit,w,&result,&error);
		printf("%10g %10g\n",x,result);
	} 
	gsl_integration_workspace_free(w);
return 0;
}
