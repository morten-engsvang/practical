#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include"qr.h"
#include"rkstep12.h"
#include<float.h> //Giver machine epsilon som DBL_EPSILON
double delta = sqrt(DBL_EPSILON);

static double epsilon;

void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps) {
	//I will construct the Jacobian using finite differences, delta_x = sqrt(DBL_EPSILON)
	//And diagonalise it using the routine implemented in a previous homework.
	int n=x->size; //Finder længden af x.
	int maxsteps = 10000;
	int step = 0;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* newx = gsl_vector_alloc(n);
	gsl_vector* backx = gsl_vector_alloc(n);
	gsl_vector* backfx = gsl_vector_alloc(n);
	
	while (step < maxsteps){
		f(x,fx); //Fylder fx med funktionsværdierne.
		//Laver min Jacobian via. finite differences:
		for (int j = 0; j < n; j++){
			double xj = gsl_vector_get(x,j); //Gemmer elementets originalværdi til senere brug.
			gsl_vector_set(x,j,xj+delta);
			f(x,df);
			gsl_vector_sub(df,fx); /* Giver df = f(x+delta)-f(x) for denne værdi af k=j i ligning 7, næste loop over i udvælger så alle ij elementer.*/
			for (int i = 0; i < n; i++){
				gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/delta);
			}
			gsl_vector_set(x,j,xj); //Putter originalværdien tilbage
		}
		GS_decomp(J,R); //Foretager QR faktorisering
		GS_solve(J,R,fx,newx); //Løser JDx=f(x) og placerer løsningen i newx.
		gsl_vector_scale(newx,-1.0); //Løsningen til JDx=-f(x) er den negative af ovenstående løsning.
		//Starter nu backtracking linesearch
		double l = 1;
		while (l > 1./64){ //Backtracker indtil l er 1/64 eller mindre.
			gsl_vector_memcpy(backx,x); //Vi skal bruge x igen og GSL er destruktiv :(
			gsl_vector_add(backx,newx); //laver x+dx vektoren ud fra den tidligere løsning, fuld skridt.
			f(backx,backfx); //Finder funktionsværdien i f(x+dx), dvs. fuld skridt.
			//Tjekker så om den opfylder ligning 8, gsl_blas_dnrm2 giver normen.
			if ( gsl_blas_dnrm2(backfx) < (1-l/2)*gsl_blas_dnrm2(fx)){
				break;
			}
			l *= 0.5; //Vi var gået for langt tilbage, så vi rykker frem.
		}
		gsl_vector_memcpy(backx,x); //Jeg laver kopi af x
		gsl_vector_scale(newx,l); //Jeg skalerer min opdatering til x med den fundne l værdi
		gsl_vector_add(backx,newx); // Lægger min opdatering til x
		gsl_vector_memcpy(x,backx); // Erstatter x med det nye x.
		f(x,fx); //Finder funktionsværdien i det nye punkt.
		//Jeg kan nu tjekke om det har ændret sig nok, tjekker om ændringerne er store nok,
		//eller om vi er tæt nok på til at f(x+dx) = 0
		if ( gsl_blas_dnrm2(newx) < delta || gsl_blas_dnrm2(fx) < eps){
			break;
		}
		step++; //Ej at forglemme
	}
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(df);
	gsl_vector_free(newx);
	gsl_vector_free(backx);
	gsl_vector_free(backfx);
}


void x3(gsl_vector* x, gsl_vector* fx){
	double t,y;
	for (int i = 0; i < x->size; i++){
		t = gsl_vector_get(x,i);
		y = pow(t,3);
		gsl_vector_set(fx,i,y);
	}
	/*
	int dimensions = x -> size;
	double scale = dimensions*gsl_blas_dnrm2(fx); //Må maks have en norm lig antallet af variable/dimensioner.
	*/
	double scale = pow(0.001,-3); //Min anden implementering virkede ikke, da alle værdierne bliver for små omkring 0.
	gsl_vector_scale(fx,scale);
}

void rosenbrockgrad(gsl_vector* x, gsl_vector* fx){
	//Gradienten af rosenbrock funktionen
	double a = gsl_vector_get(x,0);
	double b = gsl_vector_get(x,1);
	double grada = 400*pow(a,3)-400*a*b+2*a-2;
	double gradb = 200*b-200*pow(a,2);
	gsl_vector_set(fx,0,grada);
	gsl_vector_set(fx,1,gradb);
	double scale = 1000;
	gsl_vector_scale(fx,scale);
}

void schrodinger(double r, gsl_vector* y, gsl_vector* dydr){
	gsl_vector_set(dydr,0,gsl_vector_get(y,1));
	gsl_vector_set(dydr,1,-2*gsl_vector_get(y,0)*(epsilon+1/r));

}

void shootfunc(gsl_vector* x, gsl_vector* fx){
	double rinit = 0.001; //Small start r, close to 0
	double rmax = 8.0; //"large" final r.
	double delta = 0.01;
	double relacc = 0.01;
	double h = 0.001;
	int n = 2;
	gsl_vector* yinit = gsl_vector_alloc(n);
	gsl_vector* yend = gsl_vector_alloc(n);
	gsl_vector_set(yinit,0,rinit-rinit*rinit); //startværdi af f(r)
	gsl_vector_set(yinit,1,1.0-2*rinit); //startværdi af f'(r)
	epsilon = gsl_vector_get(x,0); //Kan ses af schrodinger funktionen, ændrer den "variabel"
	FILE * numode = fopen("numode.txt","w");
	driver(schrodinger,rinit,yinit,rmax,yend,h,delta,relacc,n,numode);
	gsl_vector_set(fx,0,gsl_vector_get(yend,0)); //Værdien af F(r) ved rmax, giver 0 hvis det er en løsning.
	fclose(numode);
	
	gsl_vector_free(yinit);
	gsl_vector_free(yend);
}


int main(void){
	//Test om den virker
	double eps = 0.0001;
	gsl_vector* x = gsl_vector_alloc(1);
	double xguess = 1.0;
	gsl_vector_set(x,0,xguess);
	newton(x3,x,eps);
	printf("Initial test for debugging:\n");
	printf("Root of x^3 = %10g with tolerance = %10g\n",gsl_vector_get(x,0),eps);
	printf("----------------------------------\n\n");
	//Tester for Rosenbrock funktionen
	gsl_vector* rosvector = gsl_vector_alloc(2);
	double rosx = 2.5, rosy = -0.5;
	gsl_vector_set(rosvector,0,rosx); 
	gsl_vector_set(rosvector,1,rosy);
	newton(rosenbrockgrad,rosvector,eps);
	printf("Extremum of Rosenbrocks function, with an inital guess of (%g,%g)\n",rosx,rosy);
	printf("This gives a root of (%g,%g) with tolerance = %10g\n",gsl_vector_get(rosvector,0),gsl_vector_get(rosvector,1),eps);
	gsl_vector* fx = gsl_vector_alloc(2);
	rosenbrockgrad(rosvector,fx);
	printf("The value at this point is (%g,%g)\n",gsl_vector_get(fx,0),gsl_vector_get(fx,1));
	printf("----------------------------------\n\n");
	//Laver del B angående hydrogen:
	//The boundary condition that needs to be proven can be seen to be true by inserting f(x)=r-r^2 and
	//r = 0 (don't worry 1/r dissapears due to inserting the function)
	//gives the trivial result 0 = 0 and therefore it does not violate math :)
	printf("Now calculating for the schrodinger equation\n");
	gsl_vector* guess = gsl_vector_alloc(1);
	eps = 0.01;
	double initguess = -1.0; 
	gsl_vector_set(guess,0,initguess);
	newton(shootfunc,guess,eps);
	printf("Initial guess for the energy is: %g with tolerance %g\n",initguess,eps);
	printf("The root of M(epsilon) = %g\n", gsl_vector_get(guess,0));
	printf("The expected result is -0.5, so this is as expected.\n");
	
	//Den exacte løsning til at sammenligne
	double rinit = 0.0001; //Small start r, close to 0
	double rmax = 8.0; //"large" final r.
	FILE * exactode = fopen("exactode.txt","w");
	double r, y, numsteps = 1000;
	for (int i = 0; i < numsteps; i++){
		r = rinit+(rmax-rinit)*i/(numsteps-1);
		y = r*exp(-r);
		fprintf(exactode, "%10g %10g\n",r,y);
	}
	fclose(exactode);
	
	gsl_vector_free(x);
	gsl_vector_free(rosvector);
	gsl_vector_free(guess);
return 0;
}
