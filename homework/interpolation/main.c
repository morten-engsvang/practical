#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<assert.h>
#include<stdlib.h>
#include<math.h>
#include"binsearch.h"
#include"spline.h"
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>


//Stjæler vector_print skamløst
void vector_print(char s[], gsl_vector* v){
        printf("%s\n",s);
        for (int i=0; i < v->size;i++) {
                printf("%10g ",gsl_vector_get(v,i));
        }
        printf("\n");
}



int main(int argc, char** argv){
	//Man skal give en tabel af x,y værdier og antallet af linjer som input, og værdien z,
	//der skal interpoleres
	//Jeg laver først min x og y vektorer.
	int n = atoi(argv[2]);
	double z = atof(argv[3]);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	FILE* input = fopen(argv[1],"r");
	int i = 0;
	int items;
	double tmp1;
	double tmp2;
	do{
		items = fscanf(input,"%lf %lf", &tmp1, &tmp2);
		gsl_vector_set(x,i,tmp1);
		gsl_vector_set(y,i,tmp2);
		i++;
	}
	while(i<n);
	fclose(input);
	
	FILE* output1=fopen("out_data.txt","w");
	FILE* output2=fopen("out.txt","w");	
	for(i=0;i<n-1;i++){
		double xmin=gsl_vector_get(x,i),xmax=gsl_vector_get(x,i+1);
		for(double temp=xmin;temp<=xmax;temp+=1.0/6) {
			fprintf(output1,"%10g %10g\n",temp,linterp(x,y,temp,n));
		}
	}
	double integral = linterp_integ(x,y,z,n);
	fprintf(output2,"Følgende er for lineær interpolation.\n");
	fprintf(output2,"For z=%g er det interpolerede punkt via. min implementering: %g\n",z,linterp(x,y,z,n));
	fprintf(output2,"For z=%g giver integralet fra starten af datasættet til z: %g\n",z,integral);
	fprintf(output2,"Interpoleringen kan ses i plot.svg, integralet i int_afledt.svg\n");
	
	
	double xs[n],ys[n];
	for (int i = 0; i < n; i++){
		xs[i] = gsl_vector_get(x,i);
		ys[i] = gsl_vector_get(y,i);
	}
	
	FILE* outlspline = fopen("outlspline.txt","w");
	FILE* outinteg = fopen("outlinteg.txt","w");
	gsl_interp * linear = gsl_interp_alloc(gsl_interp_linear,n);
	gsl_interp_init(linear,xs,ys,n);
	/*
	int xmin = 1;
	int xmax = 11;
	double nxmin = xmin*0.9;
	double nxmax = xmax*0.9;
	*/
	for (int i = 0; i<(2*n); i++){
		double z = 1.0+(11.0-1.0)*i/(2.0*n-1.0);
		double interp_l = gsl_interp_eval(linear,xs,ys,z,NULL);
		double integ_l=gsl_interp_eval_integ(linear,xs,ys,gsl_vector_get(x,0),z,NULL);
		fprintf(outlspline,"%10g %10g\n",z,interp_l);
		fprintf(outinteg,"%10g %10g\n",z,integ_l);
	}
	double integ_z=gsl_interp_eval_integ(linear,xs,ys,gsl_vector_get(x,0),4,NULL);
	fclose(outlspline);
	fclose(outinteg);
	fprintf(output2,"For z = %g, så giver GSL integralet som: %g\n",z,integ_z);
	
	FILE* outquadspline = fopen("outquadspline.txt","w");
	qspline * quad = quad_alloc(n, xs, ys);
	for(i=0;i<n-1;i++){
		double xmin=gsl_vector_get(x,i),xmax=gsl_vector_get(x,i+1);
		for(double temp=xmin;temp<=xmax;temp+=1.0/10) {
			fprintf(outquadspline,"%10g %10g\n",temp,quad_interp(quad,temp));
		}
	}
	
	fclose(outquadspline);
	fprintf(output2,"\nFølgende er for quadritic interpolation:\n");
	fprintf(output2,"Interpoleringen kan ses i quad.svg");
	double quad_int = quad_integ(quad,z);
	double quad_der = quad_deriv(quad,z);
	fprintf(output2,"Integralet til og med %g er %g i følge quadratic interpolation.\n",z,quad_int);
	fprintf(output2,"Den afledte i punktet er: %g",quad_der);
	fprintf(output2,"Integralet og den afledte kan findes i int_afledt.svg\n");
	
	FILE* quadinteg = fopen("quadinteg.txt","w");
	FILE* quadderiv = fopen("quadderiv.txt","w");
	FILE* lininteg = fopen("lininteg.txt","w");
	
	for (int i = 0; i<(4*n); i++){
		double z = 1.0+(11.0-1.0)*i/(4.0*n-1.0);
		double quad_integral = quad_integ(quad,z);
		double quad_derivative = quad_deriv(quad,z);
		double lin_integral = linterp_integ(x,y,z,n);
		fprintf(quadinteg,"%10g %10g\n",z,quad_integral);
		fprintf(quadderiv,"%10g %10g\n",z,quad_derivative);
		fprintf(lininteg,"%10g %10g\n",z,lin_integral);
	}
	
	fclose(quadinteg);
	fclose(quadderiv);
	fclose(lininteg);
	gsl_vector_free(x);
	gsl_vector_free(y);
	quad_free(quad);
return 0;
}
