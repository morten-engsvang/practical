#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pthread.h>
//Implementering af MC ved brug af pthreads

struct params {int n; unsigned int seed; double sum;};

void* montecarlo(void* arg){
	struct params *p = (struct params*)arg;
	p->sum=0;
	for (int i=0; i < p->n; i++) {
		double x=(double)rand_r(&(p->seed))/RAND_MAX;
		double y=(double)rand_r(&(p->seed))/RAND_MAX;
		//double z=(double)rand_r(&(p->seed))/RAND_MAX;
		//printf("Coordinates: %g, %g --> Distance from center: %g\n",x,y,z);
		if (sqrt(pow(x,2)+pow(y,2)) <= 1) {
			p->sum++;
		}
	}
	return NULL;
}

int main(void) {
	int N = 400000;
	pthread_t t1,t2,t3;
	struct params p1 = {.n=N,.seed=1,.sum=0};
	struct params p2 = {.n=N,.seed=2,.sum=0};
	struct params p3 = {.n=N,.seed=3,.sum=0};
	pthread_create(&t1,NULL,montecarlo,(void*)&p1);
	pthread_create(&t2,NULL,montecarlo,(void*)&p2);
	pthread_create(&t3,NULL,montecarlo,(void*)&p3);
	pthread_join(t1,NULL);
	pthread_join(t2,NULL);
	pthread_join(t3,NULL);
	double sum=p1.sum+p2.sum+p3.sum;
	double total=p1.n+p2.n+p3.n;
	double pi_approx = 4*sum/(total);
	printf("My approximation of pi using pthreads is: %g\nBecause the sum of points inside is: %g and total amount of points is: %g\n ",pi_approx,sum,total);
return 0;
}
