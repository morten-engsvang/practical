#include<gsl/gsl_vector.h>
#include"binsearch.h"
#include<math.h>

double linterp(gsl_vector* x, gsl_vector* y, double z, int n){
	//x,y værdierne, z er værdien jeg vil interpolere, n er længden af vektor x,y
	//Jeg identificerer starten af intervallet
	int i = binsearch(n,x,z);	
	double slope = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i)); //Jeg finder hældningen ved punktet
	return (gsl_vector_get(y,i)+slope*(z-gsl_vector_get(x,i)));//Returnerer det interpolerede punkt
}

double linear_integral(double x1, double x2, double y1, double slope){
	//Giver det lineære integral mellem to punkter givet hældningen der beskriver linjen.
	return y1*(x2-x1)+slope*pow(x2-x1,2)/2;

	//return (slope*(pow(x2,2)-pow(x1,2))*0.5+y1*(x2-x1));
}

double linterp_integ(gsl_vector* x, gsl_vector* y, double z, int n){
	//x,y værdierne, z er værdien jeg vil interpolere, n er længden af vektor x,y
	//Beregner integralet fra starten af de kendte punkter t.o.m. z
	//beregner først integralet op til punktet inden z.
	//Det gøres over alle intervallerne
	int last_known_index = binsearch(n,x,z);
	int i;
	double integral_sum = 0;
	double slope;
	for (i=0;i<last_known_index;i++) {
		slope = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
		integral_sum += linear_integral(gsl_vector_get(x,i),gsl_vector_get(x,i+1),gsl_vector_get(y,i),slope);
		//printf("Running tally for integral: %lf\n",integral_sum);
	}
	//Jeg finder så integralet fra det sidste kendte punkt ud til den interpolerede værdi.
	if (z > gsl_vector_get(x,last_known_index)){
		//printf("Z var større, da sidste kendte var: %lf\n",gsl_vector_get(x,last_known_index));
		i = last_known_index;
		slope = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
		integral_sum += linear_integral(gsl_vector_get(x,i),z,gsl_vector_get(y,i),slope);
		//printf("Running tally for integral: %lf\n",integral_sum);
	}
	return integral_sum;
}

//Definerer først en struct ligesom Fedorov
typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline * quad_alloc(int n, double *x, double *y){
	qspline * s = (qspline*) malloc(sizeof(qspline));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->b = (double*) malloc((n-1)*sizeof(double));
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc(n*sizeof(double));
	s->n = n;
	
	//Fylder min struktur:
	for (int i=0;i<n;i++){
		s->x[i]=x[i],s->y[i]=y[i];
	}
	
	//Beregner dx, ændringerne mellem x-værdierne
	//og beregner p, som er hældningen i et bestemt punkt:
	double dx[n-1], p[n-1];
	for (int i=0;i<n-1;i++) {
		dx[i]=x[i+1]-x[i]; 
		p[i]=(y[i+1]-y[i])/dx[i];
	}
	//Beregner mine c koefficienter.
	// Forlæns og baglæns rekursion som i ligning 13 og 14, hvor jeg gør det han foreslår:
	//Kører først forlæns rekursion med c[0]=0, hvor efter jeg laver baglæns rekursion fra 0.5*c[n-2]
	s->c[0]=0; //Op
	for (int i=0;i<n-2;i++) {
		s->c[i+1] = 1/dx[i+1] * (p[i+1]-p[i]-s->c[i]*dx[i]);
	}
	s->c[n-2]/=2; //Ned
	for (int i=n-3;i>=0;i--) {
		s->c[i] = 1/dx[i] * (p[i+1]-p[i]- s->c[i+1]*dx[i+1]);
	}
	//Beregner så mine b værdier: p[i]-c[i]*dx[i], defineret i ligning 15
	for(int i=0;i<n-1;i++) {
		s->b[i]=p[i]-s->c[i]*dx[i];
	}
	return s;
}

void quad_free(qspline *s){
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}

double quad_interp(qspline *s, double z){
	gsl_vector* x_vec = gsl_vector_alloc(s->n);
	for (int i = 0; i < s->n; i++){
		gsl_vector_set(x_vec,i,s->x[i]);
	}
	int i = binsearch(s->n,x_vec,z);
	double h = z-s->x[i];
	
	gsl_vector_free(x_vec);
	return s->y[i]+h*(s->b[i]+h*s->c[i]); //Giver mig det interpolerede punkt.
}

double quad_deriv(qspline * s, double z){
	gsl_vector* x_vec = gsl_vector_alloc(s->n);
	for (int i = 0; i < s->n; i++){
		gsl_vector_set(x_vec,i,s->x[i]);
	}
	int i = binsearch(s->n,x_vec,z);
	double dz = z-s->x[i];
	
	gsl_vector_free(x_vec);
	return s->b[i]+2*s->c[i]*dz; //Afledte af ligning 15. b[i]+2*c[i]*dz
}

double quad_integ(qspline * s, double z){
	gsl_vector* x_vec = gsl_vector_alloc(s->n);
	for (int i = 0; i < s->n; i++){
		gsl_vector_set(x_vec,i,s->x[i]);
	}
	int j = binsearch(s->n,x_vec,z);
	double sum = 0;
	double dx;
	//Jeg integrerer nu op til j:
	for (int i = 0; i<j; i++){
		dx = s->x[i+1]-s->x[i];
		sum += s->y[i]*dx + s->b[i]*1.0/2*pow(dx,2) + s->c[i]*1.0/3*pow(dx,3); 
		//Anti afledte af ligning 15. y[i]*dx + 1/2*b[i]*dx^2 + 1/3*c[i]*dx^3
	}
	//Vi mangler lige lidt af integralet:
	double dz = z - s->x[j];
	sum += s->y[j]*dz + s->b[j]*1.0/2*pow(dz,2) + s->c[j]*1.0/3*pow(dz,3);
	gsl_vector_free(x_vec);
	return sum; 
}
