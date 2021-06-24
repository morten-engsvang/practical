#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<float.h> //Giver machine epsilon som DBL_EPSILON
double delta = sqrt(DBL_EPSILON);

static int steplimit = 100000;

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

void higgs(gsl_vector* x, gsl_vector* fx){
	int m = 30, n = 3;
	double breitwigner;
	double value = 0;
	FILE * data = fopen("higgs.txt","r");
	gsl_matrix* datamat =	gsl_matrix_alloc(m,n);
	gsl_matrix_fscanf(data,datamat);
	double mass = gsl_vector_get(x,0);
	double gamma = gsl_vector_get(x,1);
	double A = gsl_vector_get(x,2);
	double E, sigma_i, dsigma_i;
	for (int i = 0; i < m; i++) {
		E = gsl_matrix_get(datamat,i,0);
		sigma_i = gsl_matrix_get(datamat,i,1);
		dsigma_i = gsl_matrix_get(datamat,i,2);
		breitwigner = A/(pow(E-mass,2)+pow(gamma,2)/4); 
		value += pow(breitwigner-sigma_i,2)/pow(dsigma_i,2); //Finder summen af error funktionen
	}
	gsl_vector_set(fx,0,value);
	fclose(data);
	gsl_matrix_free(datamat);
}

void plot_higgs(gsl_vector* x){
	//Thanks to Kasper Larsen for setting this up :), I could have done it myself but was too lazy.
	FILE * data = fopen("higgs.txt","r");
	FILE * plotfile = fopen("plotdata.txt","w");
	gsl_matrix * datmat = gsl_matrix_alloc(30,3);
	gsl_matrix_fscanf(data,datmat);
	
	double mass = gsl_vector_get(x,0);
	double Gamma = gsl_vector_get(x,1);
	double A = gsl_vector_get(x,2);
	double E, sigmai, dsigmai;
	
	fprintf(plotfile,"#index 0\n");
	
	for (int i = 0; i < 30; i++) {
		E = gsl_matrix_get(datmat,i,0);
		sigmai = gsl_matrix_get(datmat,i,1);
		dsigmai = gsl_matrix_get(datmat,i,2);
		fprintf(plotfile,"%10g %10g %10g\n",E,sigmai,dsigmai);
	}
	
	fprintf(plotfile,"\n\n\n#index 1\n");
	double Ein = gsl_matrix_get(datmat,0,0);
	double Efin = gsl_matrix_get(datmat,29,0);
	int numpoints = 300;
	double h = (Efin - Ein) / (numpoints - 1.0);
	for (int i = 0; i < numpoints; i++) {
		fprintf(plotfile,"%10g %10g\n",Ein+h*i,A/(pow((Ein+h*i)-mass,2)+pow(Gamma,2)/4));
	}
	gsl_matrix_free(datmat);
	fclose(plotfile);
	fclose(data);
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
	//Minimerer så for Higgs
	
	gsl_vector* higgs1 = gsl_vector_alloc(3);
	gsl_vector* higgs2 = gsl_vector_alloc(3);
	double mass = 125.3; 
	double width = 2.03; 
	double scale = 2;
	gsl_vector_set(higgs1,0,mass);
	gsl_vector_set(higgs1,1,width);
	gsl_vector_set(higgs1,2,scale);
	gsl_vector_memcpy(higgs2,higgs1);
	steps = qnewton(higgs,higgs1,eps);
	printf("Fitting the Breit-Wigner deviation function to the data gives a mass of %g, width of %g and scale-factor of %g\n",gsl_vector_get(higgs1,0),fabs(gsl_vector_get(higgs1,1)),gsl_vector_get(higgs1,2));
	printf("Starting from a value of %g, %g and %g respectively.\n",gsl_vector_get(higgs2,0),gsl_vector_get(higgs2,1),gsl_vector_get(higgs2,2));
	printf("From here it took %i steps with an accuracy of %g and a maximum step limit of %i\n",steps,eps,steplimit);
	printf("It is however very finicky");
	printf("--------------------------------------------\n\n");
	plot_higgs(higgs1);
	gsl_vector_free(x);
	gsl_vector_free(x1);
	gsl_vector_free(higgs1);
	gsl_vector_free(higgs2);
	gsl_vector_free(fx);
return 0;
}
