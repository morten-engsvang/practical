#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#define RND (double)rand()/RAND_MAX
//Den sidste definition skulle gerne give et tilfældigt tal mellem 0 og 1

void plainmc(int dim, double f(int dim,gsl_vector* x), gsl_vector* a, gsl_vector* b, int N, gsl_vector* result_error){
	//Beregner først det samlede volumen
	double V = 1.0;
	for (int i = 0; i < dim; i++){
		V *= gsl_vector_get(b,i)-gsl_vector_get(a,i);
	}
	//Gør klar til at beregne første punkt for hver dimension.
	gsl_vector* x = gsl_vector_alloc(dim);
	double sum = 0.0; //Sum af punkternes værdier
	double sum_sq = 0.0; //Sum af værdierne i anden.
	//Beregner N punkter
	for (int i = 0; i < N; i++){
		//Summerer over dimensionerne, da jeg betragter hver dimension for sig.
		//Beregner værdien i en dimension, som jeg får brug for det.
		double point = RND;
		for (int i = 0; i < dim; i++){
			gsl_vector_set(x,i, gsl_vector_get(a,i)+point*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
			double fx = f(i,x); //Funktionsværdien af punktet for i'te dimension
			sum += fx;
			sum_sq += fx*fx;
		}
	}
	double mean = sum/N; //Forventningsværdien er summen af alle funktionsværdier divideret med antallet af punkter.
	double sigma = sqrt(sum_sq/N-mean*mean); //Spredningen som defineret i ligning 4.
	gsl_vector_set(result_error,0,mean*V); //Resultat defineret i ligning 2.
	gsl_vector_set(result_error,1,sigma*V/sqrt(N)); //Fejlen som defineret i ligning 3.
	gsl_vector_free(x);
}
