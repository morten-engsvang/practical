#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include"myerf.h"


int main() {
	double xmin=-2,xmax=2;
	for(double x=xmin;x<=xmax;x+=1.0/8) {
		printf("%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),myerf(x)); 
		//Laver tabel med de værdier vi gerne vil plotte
	}
return 0;
}
