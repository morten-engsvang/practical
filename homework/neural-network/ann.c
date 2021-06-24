#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include"qnewton.h"
#include"ann.h"
static int paramcount = 3;


ann* ann_alloc(int n, double (*f)(double x), double (*f_dif)(double x), double (*f_int)(double x)){
	ann* network = malloc(sizeof(ann));
	network->n=n; //Fastsætter antallet af neuroner
	network->f=f; //Fastsætter netværkets funktion
	network->f_dif=f_dif; //Den afledte
	network->f_int=f_int; //Den integrerede, antager at hvis f passer, så passer f_dif og f_int også.
	network->params=gsl_vector_alloc(paramcount*n); //Allokerer 3 parametere pr. neuron
	return network;
}

void ann_free(ann* network){
	gsl_vector_free(network->params);
	free(network);
}

double ann_response(ann* network,double x){
	double s = 0;
	//Parameterne er allokeret som a,b,w for neuron 1 efterfulgt af a,b,w for neuron 2.
	for (int i = 0; i < network->n; i++){
		double a = gsl_vector_get(network->params, paramcount*i+0);
		double b = gsl_vector_get(network->params, paramcount*i+1);
		double w = gsl_vector_get(network->params, paramcount*i+2);
		s += network->f((x-a)/b)*w;
	}
	return s;
}

double ann_response_diff(ann* network,double x){
	double s = 0;
	for (int i = 0; i < network->n; i++){
		double a = gsl_vector_get(network->params, paramcount*i+0);
		double b = gsl_vector_get(network->params, paramcount*i+1);
		double w = gsl_vector_get(network->params, paramcount*i+2);
		s += network->f_dif((x-a)/b)*w;
	}
	return s;	
}


double ann_response_int(ann* network,double x){
	double s = 0;
	for (int i = 0; i < network->n; i++){
		double a = gsl_vector_get(network->params, paramcount*i+0);
		double b = gsl_vector_get(network->params, paramcount*i+1);
		double w = gsl_vector_get(network->params, paramcount*i+2);
		s += network->f_int((x-a)/b)*w;
	}
	return s;
}


void ann_train(ann* network, gsl_vector* xs, gsl_vector* ys){
	void cost_function(gsl_vector* p,gsl_vector* fx){
		gsl_vector_memcpy(network->params,p);
		double sum = 0;
		for (int i = 0; i < xs->size; i++){
			double xi = gsl_vector_get(xs,i);
			double yi = gsl_vector_get(ys,i);
			double fi = ann_response(network,xi);
			sum += fabs(fi-yi);
		}
		gsl_vector_set(fx,0,sum/(xs->size));
	}
	gsl_vector* p = gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);
	int steps = qnewton(cost_function,p,1e-4);
	printf("Number of steps used in optimisation used in the training was %i\n",steps);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}

