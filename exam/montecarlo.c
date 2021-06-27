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
		for (int i = 0; i < dim; i++){
			gsl_vector_set(x,i, gsl_vector_get(a,i)+RND*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
		}
		double fx = f(dim,x);
		sum += fx;
		sum_sq += fx*fx;
		
		/*
		for (int i = 0; i < dim; i++){
			gsl_vector_set(x,i, gsl_vector_get(a,i)+point*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
			double fx = f(i,x) //Funktionsværdien af punktet for i'te dimension
			sum += fx;
			sum_sq += fx*fx;
		}
		*/
	}
	double mean = sum/N; //Forventningsværdien er summen af alle funktionsværdier divideret med antallet af punkter.
	double sigma = sqrt(sum_sq/N-mean*mean); //Spredningen som defineret i ligning 4.
	gsl_vector_set(result_error,0,mean*V); //Resultat defineret i ligning 2.
	gsl_vector_set(result_error,1,sigma*V/sqrt(N)); //Fejlen som defineret i ligning 3.
	gsl_vector_free(x);
}

double vandercorput (int n, int base){
	//Implementering af van der Corput sekvensen:
	//Hvor n er 
	double q = 0;
	double bk = 1.0/base;
	while (n > 0) {
		q += (n % base)*bk;
		n /= base;
		bk /= base;
	}
	return q;
}

void halton (int n, int dim, gsl_vector* a, gsl_vector* b, gsl_vector* x){
	//Implementeringen af Halton sekvensen, generaliseringen til d dimensioner
	int base [] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,73,83,89,101};
	for (int i = 0; i < dim; i++){
		gsl_vector_set(x,i,gsl_vector_get(a,i)+vandercorput(n+1,base[i])*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
	}
	
	//gsl_vector_fprintf(stdout,x,"%g");
}

void halton2 (int n, int dim, gsl_vector* a, gsl_vector* b, gsl_vector* x){
	//Implementeringen af Halton sekvensen, generaliseringen til d dimensioner
	//int base [] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83};
	int base [] = {7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,101};
	for (int i = 0; i < dim; i++){
		gsl_vector_set(x,i,gsl_vector_get(a,i)+vandercorput(n+1,base[i])*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
	}
}

void haltonmc(int dim, double f(int dim,gsl_vector* x), gsl_vector* a, gsl_vector* b, int N, gsl_vector* result_error){
	//Beregner først det samlede volumen
	double V = 1.0;
	for (int i = 0; i < dim; i++){
		V *= gsl_vector_get(b,i)-gsl_vector_get(a,i);
	}
	//Gør klar til at beregne første punkt for hver dimension.
	gsl_vector* x = gsl_vector_alloc(dim);
	double sum1 = 0.0; //Sum af punkternes værdier
	//Beregner N punkter
	for (int i = 0; i < N/2; i++){
		//Summerer over dimensionerne, da jeg betragter hver dimension for sig.
		//Beregner værdien i en dimension, som jeg får brug for det.
		halton(i,dim,a,b,x);
		double fx = f(dim,x);
		sum1 += fx;
	}
	double sum2 = 0.0; //Sum af punkternes værdier
	for (int i = 0; i < N/2; i++){
		//Summerer over dimensionerne, da jeg betragter hver dimension for sig.
		//Beregner værdien i en dimension, som jeg får brug for det.
		halton2(i,dim,a,b,x);
		double fx = f(dim,x);
		sum2 += fx;
	}
	double mean = (sum1+sum2)/N; //Forventningsværdien er summen af alle funktionsværdier divideret med antallet af punkter.
	gsl_vector_set(result_error,0,mean*V); //Resultat defineret i ligning 2.
	gsl_vector_set(result_error,1,fabs(sum1-sum2)*V/N);//Jeg definerer fejlen på samme måde som fedorov
	gsl_vector_free(x);
}

