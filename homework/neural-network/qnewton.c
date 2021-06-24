#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<float.h> //Giver machine epsilon som DBL_EPSILON
double delta = sqrt(DBL_EPSILON);

static int steplimit = 10000;

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
