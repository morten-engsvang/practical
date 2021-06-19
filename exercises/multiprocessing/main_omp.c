#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>
#include<time.h>
//Implementering af MC ved brug af omp

int main(void) {
	double sum;
	double total;
	double local_sum;
	int num_threads;
	omp_set_num_threads(12);
	#pragma omp parallel private(local_sum) shared(sum,total,num_threads)
	{
		unsigned int seed=time(NULL);
		num_threads = 0;
		local_sum = 0;
		sum = 0;
		total = 0;
		double N = 400000;
		for (int i=0; i < N; i++) {
			double x=(double)rand_r(&seed)/RAND_MAX;
			double y=(double)rand_r(&seed)/RAND_MAX;
			if (sqrt(pow(x,2)+pow(y,2)) <= 1) {local_sum++;}
		}

		//Jeg laver et thread safe område hvor jeg joiner mine threads
		#pragma omp critical
		{
			sum += local_sum;
			total += N;
			num_threads++;

		}
		
	}
	double pi_approx = 4*sum/(total);
	printf("My approximation of pi using OMP is: %g\nBecause the sum of points inside is: %g and total amount of points is: %g\n",pi_approx,sum,total);
	printf("I used %d threads\n",num_threads);
return 0;
}
