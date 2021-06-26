#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_blas.h>
#include"decomp.h"


double functions(int i, double x){
	return pow(x,i);
	
}

void least_squares(int m, double f(int i,double x), gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c,gsl_matrix* S){
	int n = x->size;
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	gsl_matrix* R_inv = gsl_matrix_alloc(m,m);
	gsl_matrix* I = gsl_matrix_alloc(m,m);
	//Jeg laver nu matrix A og b som defineret i ligning syv i leastsq kapitlet
	for(int i=0;i<n;i++){
		double xi, yi, dyi;
		xi = gsl_vector_get(x,i);
		yi = gsl_vector_get(y,i);
		dyi = gsl_vector_get(dy,i);
		gsl_vector_set(b,i,yi/dyi);
		for(int j=0;j<m;j++){
			gsl_matrix_set(A,i,j,f(j,xi)/dyi);
		}
	}
	//Jeg skal så decompose min dannede A matrix som nu er min ønskede Q matrix
	GS_decomp(A,R);
	GS_solve(A,R,b,c);
	//Jeg har nu min løsning i form af vektor c :)
	
	//Jeg kan så også beregne covarians matricen:
	gsl_matrix_set_identity(I);
	GS_inverse(I,R,R_inv); //Udnytter min GS_inverse funktion, kan også bruges til at finde den inverse generelt :)
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,R_inv,R_inv,0,S); //Jeg har nu beregnet kovarians matricen, givet ved ligning 15.
	
	
	gsl_matrix_free(A);
	gsl_vector_free(b);
	gsl_matrix_free(R);
	gsl_matrix_free(R_inv);
	gsl_matrix_free(I);
}


int main(void){
	//Stjæler assignment af mine x og y værdier fra Fedorov:
	double x[] = {1, 2, 3, 4, 6, 9, 10, 13, 15};
	double y[] = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};
	int n=sizeof(x)/sizeof(x[0]); //Finder længden af x.
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
	gsl_matrix* S = gsl_matrix_alloc(m,m); //Kovariansmatricen 
	//Jeg kan så kalde min fit-funktion:
	least_squares(m,functions,vx,vy,vdy,c,S);
	
	gsl_vector* dc = gsl_vector_alloc(m);
	for (int k = 0; k < m; k++){
		double diag = gsl_matrix_get(S,k,k);
		gsl_vector_set(dc,k,sqrt(diag));
	}
	
	

	double c1 = gsl_vector_get(c,1);
	double dc1 = gsl_vector_get(dc,1);
	double T = -1/c1*log(2.0);
	double dT=dc1/(c1*c1); //sqrt(S_11/(c_1*c_1))
	
	double fit(double x){
		double s = 0;
		for(int k=0;k<m;k++){
			s+=gsl_vector_get(c,k)*functions(k,x);
		}
		return s;
	}
	
	double plus(int i, double x){
		return fit(x)+gsl_vector_get(dc,i)*functions(i,x);
		}

	double minus(int i, double x){
		return fit(x)-gsl_vector_get(dc,i)*functions(i,x);
		}
	
	
	//Stjæler skamløst dataudskrivning fra Fedorov:
	
	printf("Del A:\n");
	printf("Besvaret med plottet, plot.svg\n");
	printf("Jeg skulle ligeledes finde half-life af Ra-224:\n");
	printf("half-life = %.3g days\n",T);
	printf("The modern value is 3.63 days\n");
	printf("\nDel B:\n");
	printf("My value with uncertainty is: %.3g+-%.2g",T,dT);
	printf("Therefore it does not agree with the modern value.\n");
	
	FILE* data = fopen("data.txt","w");
	
	fprintf(data,"# time log(activity) delta(log(activity))\n");
	
	for(int i=0;i<n;i++)fprintf(data,"%g %g %g\n",x[i],y[i],dy[i]);
	fprintf(data,"\n\n");

	for(int i=0;i<m;i++){
		fprintf(data,"# time fit fit_plus fit_minus; k=%i\n",i);
		for(double z=x[0],dz=(x[n-1]-x[0])/64;z<=x[n-1];z+=dz)
		fprintf(data,"%g %g %g %g\n",z,fit(z),plus(i,z),minus(i,z));
	fprintf(data,"\n\n");
	}
	fclose(data);
	printf("\nDel C:\n");
	printf("Se plot.svg, har gjort samme som Fedorov, hvor jeg plotter med plus eller minus\nusikkerheden for de forskellige koefficienter\n");
	
	gsl_vector_free(vx);
	gsl_vector_free(vy);
	gsl_vector_free(vdy);
	gsl_vector_free(c);
	gsl_matrix_free(S);
	gsl_vector_alloc(m);
return 0;
}
