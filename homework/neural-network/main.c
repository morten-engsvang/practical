#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<math.h>
#include"ann.h"

double activation_function(double x){
	return x*exp(-x*x);
}

double active_dif(double x){
	return exp(-x*x)*(1-2*x*x);
}

double active_int(double x){
	return -exp(-x*x)*0.5;
}

double test_function(double x){
	return cos(5*x-1)*exp(-x*x);
}

int main(void){
	//Initialiserer et netværk af n neuroner, med den tidligere definerede activation function
	//som jeg har valgt til at være en wavelet.
	int n = 5;
	ann* network = ann_alloc(n,activation_function,active_dif,active_int);
	double a = -3.14, b = 3.14;
	int nx = 50;
	gsl_vector* vx = gsl_vector_alloc(nx);
	gsl_vector* vy = gsl_vector_alloc(nx);
	
	//Laver punkterne som den skal trænes på:
	for (int i = 0; i < nx; i++){
		double x = a+(b-a)*i/(nx-1);
		double f = test_function(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
	}
	
	//Fastsætter startværdier for parameterne
	for (int i = 0; i < network->n; i++){
		gsl_vector_set(network->params,3*i+0,a+(b-a)*i/((network->n)-1));
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,2);
	}
	//Vi kan så træne:
	ann_train(network,vx,vy);
	FILE * function = fopen("function.txt","w");
	FILE * neural = fopen("neural.txt","w");
	FILE * neural_dif = fopen("neural_dif.txt","w");
	FILE * neural_int = fopen("neural_int.txt","w");
	FILE * params = fopen("params.txt","w");
	for (int i = 0; i < vx->size; i++){
		double x=gsl_vector_get(vx,i);
		double f=gsl_vector_get(vy,i);
		fprintf(function,"%g %g\n",x,f);
	}
	double dz = 1.0/64;
	for (double z = a; z <= b; z += dz){
		double y = ann_response(network,z);
		double dif = ann_response_diff(network,z);
		double integ = ann_response_int(network,z);
		fprintf(neural,"%g %g\n",z,y);
		fprintf(neural_dif,"%g %g\n",z,dif);
		fprintf(neural_int,"%g %g\n",z,integ);
	}
	printf("Training is tested for the function f(x)=cos(5x-1)*exp(-x*x) from -1 to 1\n");
	printf("The result after training is printed in plit.svg\n");
	printf("Final parameters are printed in params.txt\n");
	for (int i = 0; i < network->n; i++){
		double ai = gsl_vector_get(network->params,3*i+0);
		double bi = gsl_vector_get(network->params,3*i+1);
		double wi = gsl_vector_get(network->params,3*i+2);
		fprintf(params,"i=%i ai,bi,wi = %g %g %g\n",i,ai,bi,wi);
	}
	
	fclose(function);
	fclose(neural);
	fclose(neural_dif);
	fclose(neural_int);
	fclose(params);
	gsl_vector_free(vx);
	gsl_vector_free(vy);
	ann_free(network);
return 0;
}
