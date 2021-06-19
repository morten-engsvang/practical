#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<assert.h>
#include<stdlib.h>
#include<math.h>
#include"binsearch.h"
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
		//printf("Hello, i = %d\n",i);
		gsl_vector_set(x,i,tmp1);
		gsl_vector_set(y,i,tmp2);
		//printf("x: %lf og y: %lf\n",tmp1,tmp2);
		i++;
	}
	while(i<n);
	fclose(input);
	/*
	printf("I'm free!\n");
	//Hurtigt debug:
	vector_print("Tjekker x:",x);
	vector_print("Tjekker y:",y);
	double z = 3;
	double tjek = linterp(x,y,z,n);
	printf("For z=%lf får jeg punktet: %lf\n",z,tjek);
	*/
	//Jeg kan derefter integrere op til punktet z.
	
	//double integral = linterp_integ(x,y,z,n);
	/*
	printf("For z=%lf får jeg integralet fra start t.o.m. z: %lf\n",z,integral);
	*/
	FILE* output1=fopen("out_data.txt","w");
	FILE* output2=fopen("out.txt","w");	
	for(i=0;i<n-1;i++){
		double xmin=gsl_vector_get(x,i),xmax=gsl_vector_get(x,i+1);
		for(double temp=xmin;temp<=xmax;temp+=1.0/12) {
			fprintf(output1,"%10g %10g\n",temp,linterp(x,y,temp,n));
		}
	}
	double integral = linterp_integ(x,y,z,n);
	fprintf(output2,"For z=%lf er det interpolerede punkt via. min implementering: %lf\n",z,linterp(x,y,z,n));
	fprintf(output2,"For z=%lf giver integralet fra starten af datasættet til z: %lf\n",z,integral);
	//gsl_interp* w;
	//w = gsl_interp_alloc(gsl_interp_linear,n);
	//double interp = gsl_interp_eval(w,x,y,z);
	//fprintf(output2,"De indbyggede funktioner i gsl giver: %lf",interp);
	gsl_vector_free(x);
	gsl_vector_free(y);
return 0;
}
