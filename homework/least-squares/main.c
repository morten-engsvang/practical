#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_blas.h>
#include"decomp.h"

double functions(int i, double x){
	//I bestemmer exponenten af mine exponentielle funktioner
	//x er punktet hvor funktionsværdien ønskes.
	double result = pow(x,i);
	return result;
	/*
	switch(i){
		case 0: 
			return 1; 
			break;
		case 1:
			return x;
			break;
		case 2:
			return x*x;
		       	break;
		case 3:
			return pow(x,3);
			break;
		case 4:
			return pow()
		default:
			return NAN;
	}
	*/
}

void least_squares(int m,gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c){
	int n = c->size;
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	gsl_matrix* R_inv = gsl_matrix_alloc(m,m);
	//Jeg laver nu matrix A og b som defineret i ligning syv i leastsq kapitlet
	for(int i=0;i<n;i++){
		double x_i, y_i, dy_i;
		x_i = gsl_vector_get(x,i);
		y_i = gsl_vector_get(y,i);
		dy_i = gsl_vector_get(dy,i);
		gsl_vector_set(b,i,y_i/dy_i);
		for(int j=0;j<m;j++){
			gsl_matrix_set(A,i,j,functions(j,x_i)/dy_i);
		}
	}
	//Jeg skal så decompose min dannede A matrix som nu er min ønskede Q matrix
	GS_decomp(A,R);
	GS_solve(A,R,b,c);
	//Jeg har nu min løsning i form af vektor c :)

	gsl_matrix_free(A);
	gsl_vector_free(b);
	gsl_matrix_free(R);
	gsl_matrix_free(R_inv);
}


int main(void){
	//Stjæler assignment af mine x og y værdier fra Fedorov:
	double x[] = {1, 2, 3, 4, 6, 9, 10, 13, 15};
	double y[] = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};
	int n = 9;
	//Finder først normal usikkerhed:
	double dy[n];
	for(int i=0;i<n;i++){
		dy[i]=0.05*y[i];
	}
	//Omdanner y til logaritme og dy til logaritmisk usikkerhed:
	for(int i=0;i<n;i++){
		dy[i]/=y[i];
		y[i]=log(y[i]);
	}
	//Opsætter mine vektorer:
	gsl_vector* vx = gsl_vector_alloc(n);
	gsl_vector* vy = gsl_vector_alloc(n);
	gsl_vector* vdy = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(vx,i,x[i]);
		gsl_vector_set(vy,i,y[i]);
		gsl_vector_set(vdy,i,dy[i]);
	}
	//m bestemmer antallet af funktioner jeg vil medtage i mit fit:
	int m = 2; //Skal kun bruge to, den konstante og så den i første.
	gsl_vector* c = gsl_vector_alloc(m);
	//gsl_matrix* S = gsl_matrix_alloc(m,m); //Kovariansmatricen hvis jeg kommer så langt.
	
	//Jeg kan så kalde min fit-funktion:
	least_squares(m,vx,vy,vdy,c);
	//printf("Jeg får følgende koefficient vektor c=\n");
	//gsl_vector_fprintf(stdout,c,"%g");
	
	double fit(double x){
		double s = 0;
		for(int k=0;k<m;k++){
			s+=gsl_vector_get(c,k)*functions(k,x);
		}
		return s;
	}

	double c1 = gsl_vector_get(c,1);
	double T = -1/c1*log(2.0);
	
	//Stjæler skamløst dataudskrivning fra Fedorov:
	printf("# half-life = %.3g days\n",T);
	
	printf("# time log(activity) delta(log(activity))\n");
	for(int i=0;i<n;i++)printf("%g %g %g\n",x[i],y[i],dy[i]);
	printf("\n\n");

	printf("# time fit\n");
	for(double z=x[0],dz=(x[n-1]-x[0])/64;z<=x[n-1];z+=dz){
		printf("%g %g\n",z,fit(z));
	}
	printf("\n\n");
return 0;
}
