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

static int steplimit = 100000;

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

void numeric_gradient(void F(gsl_vector* x,gsl_vector* fx), gsl_vector* x, gsl_vector* gradient){
	//Sætter den givne vector gradient til gradienten i punktet givet ved vektor x for 
	//funktionen F, uden at ændre på vektor x.
	int n=x->size;
	gsl_vector* fx = gsl_vector_alloc(1);
	gsl_vector* fx1 = gsl_vector_alloc(1);
	F(x,fx); //Finder funktionsværdierne i punktet
	for (int i = 0; i < n; i++){
		double dx = delta;
		double xi = gsl_vector_get(x,i);
		if (fabs(xi)>sqrt(delta)){
			dx = fabs(xi)*delta;
		}
		gsl_vector_set(x,i,xi+dx);
		F(x,fx1);
		double grad = (gsl_vector_get(fx1,0)-gsl_vector_get(fx,0))/dx; //Numerisk værdi af gradienten
		gsl_vector_set(gradient,i,grad);
		gsl_vector_set(x,i,xi); //Resetter x til udgangspunktet.
	}
}


int qnewton(void F(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps){
		//F er funktionen, der ønskes minimeret
		//X er startpunktet, output som minimumspunktet.
		//Den numeriske gradienten kommer fra funktionen ovenover numeric_gradient
		//Bruger en rank-1 update i form af den symmetriske broyden opdatering
		//Og vi har back-tracking linesearch :)
		
		int nsteps = 0;
		int n = x->size;
		gsl_matrix* B = gsl_matrix_alloc(n,n);
		gsl_vector* step = gsl_vector_alloc(n);
		gsl_vector* z = gsl_vector_alloc(n);
		gsl_vector* gradient = gsl_vector_alloc(n);
		gsl_vector* grad_z = gsl_vector_alloc(n);
		gsl_vector* y = gsl_vector_alloc(n);
		gsl_vector* u = gsl_vector_alloc(n);
		gsl_vector* a = gsl_vector_alloc(n);
		gsl_vector* fx = gsl_vector_alloc(1);
		gsl_vector* fz = gsl_vector_alloc(1);
		gsl_matrix_set_identity(B);
		numeric_gradient(F,x,gradient); //Beregner gradienten i startpunktet.
		F(x,fx); //Finder funktionsværdien i startpunktet.
		while (nsteps < steplimit){
			//Beregner normen af gradienten og placerer den i double norm2gx.
			double norm2gx = 0;
			for (int i = 0; i < n; i++){
				norm2gx += gsl_vector_get(gradient,i)*gsl_vector_get(gradient,i);
			}
			if (norm2gx < eps*eps) {
				break;
			}
			nsteps++;
			//Beregner skridet frem:
			for (int i = 0; i < n; i++){
				double step_change = 0;
				for (int j = 0; j < n; j++){
					step_change -= gsl_matrix_get(B,i,j)*gsl_vector_get(gradient,j);
				}
				gsl_vector_set(step,i,step_change);
			}
			double lambda = 1;
			while (1){
				for (int i = 0; i < n; i++){
					gsl_vector_set(z,i,gsl_vector_get(x,i)+gsl_vector_get(step,i));
				}
				F(z,fz); //Beregner funktionsværdien i det næste punkt, z.
				double sTg = 0;
				for (int i = 0; i < n; i++){
					sTg += gsl_vector_get(step,i)*gsl_vector_get(gradient,i);
				}
				if (gsl_vector_get(fz,0) < gsl_vector_get(fx,0)+0.01*sTg) {
					break;
				}
				if (lambda < delta) {
					gsl_matrix_set_identity(B);
					break;
				}
				lambda /= 2;
				for (int i = 0; i < n; i++){
					double prev_i = gsl_vector_get(step,i);
					prev_i /= 2;
					gsl_vector_set(step,i,prev_i);
				}
			}
			numeric_gradient(F,x,grad_z); //Beregner gradienten i det nye punkt.
			//Beregner y, som er ændringen i gradienten mellem punkterne x og z.
			for (int i = 0; i < n; i++){
				gsl_vector_set(y,i,gsl_vector_get(grad_z,i)-gsl_vector_get(gradient,i));
			}
			for (int i = 0; i < n; i++){
				double change = 0;
				for (int j = 0; j < n; j++){
					change -= gsl_matrix_get(B,i,j)*gsl_vector_get(y,j);
				}
				gsl_vector_set(u,i,gsl_vector_get(step,i)+change);
			}
			double sTy = 0;
			for (int i = 0; i < n; i++){
				sTy += gsl_vector_get(step,i)*gsl_vector_get(y,i);
			}
			//Laver nu den symmetriske Broyden update
			if (fabs(sTy) > delta) {
				double uTy = 0;
				for (int i = 0; i < n; i++){
					uTy += gsl_vector_get(u,i)*gsl_vector_get(y,i);
				}
				double gamma = uTy/(2*sTy);
				for (int i = 0; i < n; i++){
					gsl_vector_set(a,i,gsl_vector_get(u,i)-gamma*gsl_vector_get(step,i));
				}
				for (int i = 0; i < n; i++){
					for (int j = 0; j < n; j++){
						double B_prev = gsl_matrix_get(B,i,j);
						B_prev += gsl_vector_get(a,i)*gsl_vector_get(step,j)/sTy;
						B_prev += gsl_vector_get(step,i)*gsl_vector_get(a,j)/sTy;
						gsl_matrix_set(B,i,j,B_prev);
					}
				}
				
			}
			for (int i = 0; i < n; i++){
				gsl_vector_set(x,i,gsl_vector_get(z,i));
				gsl_vector_set(gradient,i,gsl_vector_get(grad_z,i));
			}
			gsl_vector_set(fx,0,gsl_vector_get(fz,0));
		}
		
		gsl_matrix_free(B);
		gsl_vector_free(step);
		gsl_vector_free(z);
		gsl_vector_free(gradient);
		gsl_vector_free(grad_z);
		gsl_vector_free(u);
		gsl_vector_free(y);
		gsl_vector_free(a);
		return nsteps;
}

void rosenbrock(gsl_vector* x, gsl_vector* fx){
	double a = gsl_vector_get(x,0);
	double b = gsl_vector_get(x,1);
	double value = pow(1-a,2)+100*pow(b-pow(a,2),2);
	gsl_vector_set(fx,0,value);
}

void himmelblau(gsl_vector* x, gsl_vector* fx){
	double a = gsl_vector_get(x,0);
	double b = gsl_vector_get(x,1);
	double value = pow(a*a+b-11,2)+pow(a+b*b-7,2);
	gsl_vector_set(fx,0,value);
}


int main(void){
	double eps = 0.00001;
	gsl_vector* x = gsl_vector_alloc(2);
	gsl_vector* x1 = gsl_vector_alloc(2);
	gsl_vector* fx = gsl_vector_alloc(1);
	//Beregner først for Rosenbrock:
	gsl_vector_set(x,0,1);
	gsl_vector_set(x,1,-1);
	gsl_vector_memcpy(x1,x);
	int steps = qnewton(rosenbrock,x,eps);
	printf("A minimum of the Rosenbrock's valley function is (%g,%g)\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	rosenbrock(x,fx);
	printf("The value at this point is %g\n",gsl_vector_get(fx,0));
	printf("Starting from the point (%g,%g) it took %i steps.\n",gsl_vector_get(x1,0),gsl_vector_get(x1,1),steps);
	printf("With an accuracy of %g and a maximum step limit of %i\n",eps,steplimit);
	printf("--------------------------------------------\n\n");
	//Beregner så for himmelblau:
	gsl_vector_set(x,0,1);
	gsl_vector_set(x,1,-1);
	gsl_vector_memcpy(x1,x);
	steps = qnewton(himmelblau,x,eps);
	printf("A minimum of the Himmelblau's function is (%g,%g)\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	himmelblau(x,fx);
	printf("The value at this point is %g\n",gsl_vector_get(fx,0));
	printf("Starting from the point (%g,%g) it took %i steps.\n",gsl_vector_get(x1,0),gsl_vector_get(x1,1),steps);
	printf("With an accuracy of %g and a maximum step limit of %i\n",eps,steplimit);
	printf("--------------------------------------------\n\n");
	
	
	gsl_vector_free(x);
	gsl_vector_free(x1);
	gsl_vector_free(fx);
return 0;
}
