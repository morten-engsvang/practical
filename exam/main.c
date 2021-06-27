#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include"montecarlo.h"
#include"minimum.h"

static int points = 0; //Giver mig det endelige antal af punkter

static int flag = 0; //Flag til at bestemme søjle jeg skal bruge fra plot.txt

double function(double x){
	return 1.0/sqrt(x);
}

void errfunc(gsl_vector* x, gsl_vector* fx){
	//x som input er en scaling factor til min funktion, dvs. længde 1
	//fx skal være den samlede fejl, dvs. længde 1
	FILE * data = fopen("plot.txt","r");
	gsl_matrix* datamat = gsl_matrix_alloc(points,5);
	gsl_matrix_fscanf(data,datamat);
	//printf("Vi kommer ind i errfunc\n");
	//printf("Flag = %i og antal punkter er %i\n",flag,points);
	double toterr = 0;
	double montecarloerr = 0;
	double functionvalue = 0;
	for (int i = 0; i < points; i++){
		montecarloerr = gsl_matrix_get(datamat,i,flag);
		functionvalue = gsl_vector_get(x,1)+gsl_vector_get(x,0)*function(gsl_matrix_get(datamat,i,0));
		toterr += fabs(montecarloerr-functionvalue);
	}
	//printf("Total error is right now: %g\n",toterr);
	fclose(data);
	gsl_vector_set(fx,0,toterr);
}

void errfunc2(gsl_vector* x, gsl_vector* fx){
	//x som input er en scaling factor til min funktion, dvs. længde 1
	//fx skal være den samlede fejl, dvs. længde 1
	FILE * data = fopen("plot2.txt","r");
	gsl_matrix* datamat = gsl_matrix_alloc(points,5);
	gsl_matrix_fscanf(data,datamat);
	//printf("Vi kommer ind i errfunc\n");
	//printf("Flag = %i og antal punkter er %i\n",flag,points);
	double toterr = 0;
	double montecarloerr = 0;
	double functionvalue = 0;
	for (int i = 0; i < points; i++){
		montecarloerr = gsl_matrix_get(datamat,i,flag);
		functionvalue = gsl_vector_get(x,1)+gsl_vector_get(x,0)*function(gsl_matrix_get(datamat,i,0));
		toterr += fabs(montecarloerr-functionvalue);
	}
	//printf("Total error is right now: %g\n",toterr);
	fclose(data);
	gsl_vector_set(fx,0,toterr);
}


double integral1(int dim, gsl_vector* x){
	double temp = 1-cos(gsl_vector_get(x,0))*cos(gsl_vector_get(x,1))*cos(gsl_vector_get(x,2));
	double result = pow(temp,-1.0)/pow(M_PI,3);
	return result;
}

double integral2(int dim, gsl_vector* x){
	return 4*sqrt(1-pow(gsl_vector_get(x,0),2));
}


int main(void){
	//Jeg undersøger det først for integral1
	printf("The first integral to be attempted will be the integral given in exercise A of the Monte Carlo Integration.\n");
	printf("This will be done for the Plain Monte Carlo algorithm and the Quasi-random Monte Carlo algorithm using the Halton sequence\n");
	printf("For these I will calculate the error, for a given amount of sampling points, to which I will attempt to fit the function f(x)=a+b*1/sqrt(x)\n");
	printf("This is because the error from a random walk is expected to converge with 1/sqrt(N), where N is the amount of sampling points used.\n");
	printf("I use both the estimated error and the actual error such that I have two curves pr. algorithm pr. integral.\n");
	printf("The fitting is done by defining a function that returns the error between the fitting function and the points and minimising it using\nthe quasi-newton algorithm from the Minimisation exercise.\n");
	printf("The result of this can be found in plot_plain.svg and plot_halton.svg\n");
	//Finder først fejlen som funktion af mængden af sampling points
	int dim = 3;
	gsl_vector* a = gsl_vector_alloc(dim);
	gsl_vector* b = gsl_vector_alloc(dim);
	gsl_vector* result_error = gsl_vector_alloc(2);
	for (int i = 0; i < dim; i++){
		gsl_vector_set(a,i,0);
		gsl_vector_set(b,i,M_PI);
	}
	
	double actual = 1.3932039296856768591842462603255;
	FILE * data = fopen("plot.txt","w");
	
	for (int N = 7; N < 10; N += 1){
		plainmc(dim,integral1,a,b,N,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral1,a,b,N,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data, "%i %10g %10g %10g %10g\n", N, p_est_error, p_act_error, h_est_error, h_act_error);
		points++;
	}
	
	for (int N = 10; N < 100; N += 10){
		plainmc(dim,integral1,a,b,N,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral1,a,b,N,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data, "%i %10g %10g %10g %10g\n", N, p_est_error, p_act_error, h_est_error, h_act_error);
		points++;
	}
	
	for (int N = 100; N < 1000; N += 50){
		plainmc(dim,integral1,a,b,N,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral1,a,b,N,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data, "%i %10g %10g %10g %10g\n", N, p_est_error, p_act_error, h_est_error, h_act_error);
		points++;
	}
	
	for (int N = 1000; N <= 2000; N += 100){
		plainmc(dim,integral1,a,b,N,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral1,a,b,N,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data, "%i %10g %10g %10g %10g\n", N, p_est_error, p_act_error, h_est_error, h_act_error);
		points++;
	}
	printf("\nStart of technical information:\n");
	printf("Number of datapoints: %i\n",points);
	fclose(data);
	gsl_vector* x = gsl_vector_alloc(2); //Den første er scaling, den anden er offset
	gsl_vector* scalings = gsl_vector_alloc(4);
	gsl_vector* offsets = gsl_vector_alloc(4);
	gsl_vector* stepused = gsl_vector_alloc(4);
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	//printf("The scaling value starts at %g\n",gsl_vector_get(x,0));
	double eps = 0.1;
	
	//Jeg kan nu fitte 1/sqrt(N):
	FILE * fit = fopen("fit.txt","w");
	int steps = 0;
	
	flag = 1; //Jeg bruger plainmc estimated error
	steps = qnewton(errfunc,x,eps);
	gsl_vector_set(scalings,flag-1,gsl_vector_get(x,0));
	gsl_vector_set(offsets,flag-1,gsl_vector_get(x,1));
	gsl_vector_set(stepused,flag-1,steps);
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	
	flag = 2; //Jeg bruger plainmc actual error
	steps = qnewton(errfunc,x,eps);
	gsl_vector_set(scalings,flag-1,gsl_vector_get(x,0));
	gsl_vector_set(offsets,flag-1,gsl_vector_get(x,1));
	gsl_vector_set(stepused,flag-1,steps);
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	
	flag = 3; //Jeg bruger hamiltonmc estimated error
	steps = qnewton(errfunc,x,eps);
	gsl_vector_set(scalings,flag-1,gsl_vector_get(x,0));
	gsl_vector_set(offsets,flag-1,gsl_vector_get(x,1));
	gsl_vector_set(stepused,flag-1,steps);
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	
	flag = 4; //Jeg bruger hamiltonmc actual error
	steps = qnewton(errfunc,x,eps);
	gsl_vector_set(scalings,flag-1,gsl_vector_get(x,0));
	gsl_vector_set(offsets,flag-1,gsl_vector_get(x,1));
	gsl_vector_set(stepused,flag-1,steps);
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	
	double plainestimate = 0;
	double plainactual = 0;
	double haltonestimate = 0;
	double haltonactual = 0;
	for (int N = 7; N < 2000; N += 10){
		plainestimate = fabs(gsl_vector_get(scalings,0))*function(N);
		plainactual = fabs(gsl_vector_get(scalings,1))*function(N);
		haltonestimate = fabs(gsl_vector_get(scalings,2))*function(N);
		haltonactual = fabs(gsl_vector_get(scalings,3))*function(N);
		fprintf(fit,"%i %10g %10g %10g %10g \n", N, plainestimate, plainactual, haltonestimate, haltonactual);
	}
	
	fclose(fit);
	
	printf("The scaling values are %g %g %g %g\n",fabs(gsl_vector_get(scalings,0)),fabs(gsl_vector_get(scalings,1)),fabs(gsl_vector_get(scalings,2)),fabs(gsl_vector_get(scalings,3)));
	
	printf("The offset values are %g %g %g %g\n",gsl_vector_get(offsets,0),gsl_vector_get(offsets,1),gsl_vector_get(offsets,2),gsl_vector_get(offsets,3));
	
	printf("The steps used are %g %g %g %g\n\n",gsl_vector_get(stepused,0),gsl_vector_get(stepused,1),gsl_vector_get(stepused,2),gsl_vector_get(stepused,3));
	
	printf("It can be observed that the fit is far from optimal especially for larger N, which is primarily due to the large error values for small N which can give very large errors.\n");
	
	printf("This is even the case when I start at N = 7 as I have done here.\n");
	
	printf("\nThe second integral I will be looking at is a 1D integral: 4*sqrt(1-x^2) from 0 to 1 which is supposed to give pi\n");
	
	dim = 1;
	gsl_vector* a2 = gsl_vector_alloc(dim);
	gsl_vector* b2 = gsl_vector_alloc(dim);
	gsl_vector_set(a2,0,0);
	gsl_vector_set(b2,0,1);
	actual = M_PI;
	points = 0;
	FILE * data2 = fopen("plot2.txt","w");
	
	for (int N = 3; N < 10; N += 1){
		plainmc(dim,integral2,a2,b2,N,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral2,a2,b2,N,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data2, "%i %10g %10g %10g %10g\n", N, p_est_error, p_act_error, h_est_error, h_act_error);
		points++;
	}
	
	for (int N = 10; N < 100; N += 10){
		plainmc(dim,integral2,a2,b2,N,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral2,a2,b2,N,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data2, "%i %10g %10g %10g %10g\n", N, p_est_error, p_act_error, h_est_error, h_act_error);
		points++;
	}
	
	for (int N = 100; N < 1000; N += 50){
		plainmc(dim,integral2,a2,b2,N,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral2,a2,b2,N,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data2, "%i %10g %10g %10g %10g\n", N, p_est_error, p_act_error, h_est_error, h_act_error);
		points++;
	}
	
	for (int N = 1000; N <= 2000; N += 100){
		plainmc(dim,integral2,a2,b2,N,result_error);
		double p_est_error = gsl_vector_get(result_error,1);
		double p_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		haltonmc(dim,integral2,a2,b2,N,result_error);
		double h_est_error = gsl_vector_get(result_error,1);
		double h_act_error = fabs(gsl_vector_get(result_error,0)-actual);
		fprintf(data2, "%i %10g %10g %10g %10g\n", N, p_est_error, p_act_error, h_est_error, h_act_error);
		points++;
	}
	fclose(data2);
	printf("\nStart of technical information:\n");
	printf("Number of datapoints: %i\n",points);
	
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	//printf("The scaling value starts at %g\n",gsl_vector_get(x,0));
	eps = 0.1;
	
	//Jeg kan nu fitte 1/sqrt(N):
	FILE * fit2 = fopen("fit2.txt","w");
	steps = 0;
	
	flag = 1; //Jeg bruger plainmc estimated error
	steps = qnewton(errfunc2,x,eps);
	gsl_vector_set(scalings,flag-1,gsl_vector_get(x,0));
	gsl_vector_set(offsets,flag-1,gsl_vector_get(x,1));
	gsl_vector_set(stepused,flag-1,steps);
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	
	flag = 2; //Jeg bruger plainmc actual error
	steps = qnewton(errfunc2,x,eps);
	gsl_vector_set(scalings,flag-1,gsl_vector_get(x,0));
	gsl_vector_set(offsets,flag-1,gsl_vector_get(x,1));
	gsl_vector_set(stepused,flag-1,steps);
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	
	flag = 3; //Jeg bruger hamiltonmc estimated error
	steps = qnewton(errfunc2,x,eps);
	gsl_vector_set(scalings,flag-1,gsl_vector_get(x,0));
	gsl_vector_set(offsets,flag-1,gsl_vector_get(x,1));
	gsl_vector_set(stepused,flag-1,steps);
	gsl_vector_set(x,0,20);
	gsl_vector_set(x,1,0);
	
	flag = 4; //Jeg bruger hamiltonmc actual error
	steps = qnewton(errfunc2,x,eps);
	gsl_vector_set(scalings,flag-1,gsl_vector_get(x,0));
	gsl_vector_set(offsets,flag-1,gsl_vector_get(x,1));
	gsl_vector_set(stepused,flag-1,steps);
	
	plainestimate = 0;
	plainactual = 0;
	haltonestimate = 0;
	haltonactual = 0;
	for (int N = 3; N < 2000; N += 10){
		plainestimate = fabs(gsl_vector_get(scalings,0))*function(N);
		plainactual = fabs(gsl_vector_get(scalings,1))*function(N);
		haltonestimate = fabs(gsl_vector_get(scalings,2))*function(N);
		haltonactual = fabs(gsl_vector_get(scalings,3))*function(N);
		fprintf(fit2,"%i %10g %10g %10g %10g \n", N, plainestimate, plainactual, haltonestimate, haltonactual);
	}
	
	fclose(fit2);
	
	printf("The scaling values are %g %g %g %g\n",fabs(gsl_vector_get(scalings,0)),fabs(gsl_vector_get(scalings,1)),fabs(gsl_vector_get(scalings,2)),fabs(gsl_vector_get(scalings,3)));
	
	printf("The offset values are %g %g %g %g\n",gsl_vector_get(offsets,0),gsl_vector_get(offsets,1),gsl_vector_get(offsets,2),gsl_vector_get(offsets,3));
	
	printf("The steps used are %g %g %g %g\n\n",gsl_vector_get(stepused,0),gsl_vector_get(stepused,1),gsl_vector_get(stepused,2),gsl_vector_get(stepused,3));
	
	printf("This is integral is much better behaved, and it can be seen that the fit of the 1/sqrt(x) is much better at describing the trend.\n");
	
	
	gsl_vector_free(a);
	gsl_vector_free(b);
	gsl_vector_free(result_error);
	gsl_vector_free(x);
	gsl_vector_free(scalings);
	gsl_vector_free(offsets);
	gsl_vector_free(stepused);
	gsl_vector_free(a2);
	gsl_vector_free(b2);
return 0;
}
