#include<gsl/gsl_vector.h>
#include<math.h>

void vector_print(char s[], gsl_vector* v){
        printf("%s\n",s);
        for (int i=0; i < v->size;i++) {
                printf("%10g ",gsl_vector_get(v,i));
        }
        printf("\n");
}


void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err, int n){
/* Jeg implementerer mid-point metoden. */
gsl_vector* k0 = gsl_vector_alloc(n);
gsl_vector* dydt = gsl_vector_alloc(n);
gsl_vector* y12 = gsl_vector_alloc(n);
gsl_vector* temp = gsl_vector_alloc(n);
gsl_vector* k = gsl_vector_alloc(n);
gsl_vector* err_k = gsl_vector_alloc(n);

//printf("Jeg opdaterer den afledte i det første punkt\n");

f(t,yt,dydt); /* Opdaterer dydt så det er de afledte i det første punkt */
gsl_vector_memcpy(k0,dydt); /*  Vi sætter k0 lig de afledte af funktionen i det første punkt*/
/*  Jeg gør så klar til at beregne k12  */

//printf("Jeg beregner k12\n");

double x12 = t + 0.5*h;
gsl_vector_memcpy(y12,yt);
gsl_vector_memcpy(temp,k0);
gsl_vector_scale(temp,0.5*h);
gsl_vector_add(y12,temp);

//printf("Jeg opdater den afledte i midtpunktet\n");

f(x12,y12,dydt); /* Opdaterer dydt med midtpunktsværdierne, dvs. beregner k12 */
gsl_vector_memcpy(k,dydt); /* Sætter k = k12 */
gsl_vector_memcpy(err_k,k); /* Bruges til error estimate */
/* Jeg beregner nu det næste skridt yh*/
gsl_vector_memcpy(yh,yt);
gsl_vector_scale(k,h);
gsl_vector_add(yh,k); /* Endelige bud på yh */
/* Og til sidst error estimate */

//printf("Beregner error estimate\n");

gsl_vector_add(k0,err_k);
gsl_vector_scale(k0,0.5*h);

//printf("Endelige estimat er beregnet\n");

gsl_vector_memcpy(err,k0); /* Endelige fejl estimat */


gsl_vector_free(k0);
gsl_vector_free(dydt);
gsl_vector_free(y12);
gsl_vector_free(temp);
gsl_vector_free(k);
gsl_vector_free(err_k);
}

void driver(void f(double t ,gsl_vector* y,gsl_vector* dydt), double a, gsl_vector* ya, double b, gsl_vector* yb, double h, double acc, double eps, int n, FILE * outstream) {
	if (n == 2){
		fprintf(outstream,"%g %g %g\n",a,gsl_vector_get(ya,0),gsl_vector_get(ya,1));
	}
	if (n == 3){
		fprintf(outstream,"%g %g %g %g\n",a,gsl_vector_get(ya,0),gsl_vector_get(ya,1),gsl_vector_get(ya,2));
	}
	gsl_vector* yh = gsl_vector_alloc(n);
	gsl_vector* err = gsl_vector_alloc(n);
	double h_new = h;
	int done = 0; // Flag to signal completion
	//printf("Starter driver loop\n");
	while(done == 0){
		/* Tjekker om h er for langt */
		if (h_new > b-a){
			h_new = b-a;
		}		
				
		/* Forsøger et trin frem */
		//printf("Forsøger et trin frem\n");
		//printf("---------------------------\n");
		rkstep12(f,a,ya,h_new,yh,err,n);
		//printf("---------------------------\n");
		//printf("Lykkedes med at gå et trin frem\n");
		double s = 0;
		for (int i=0;i<n;i++){
			s += pow(gsl_vector_get(err,i),2);
		}
		double ei = sqrt(s);
		s = 0;
		for (int i=0;i<n;i++){
			s += pow(gsl_vector_get(yh,i),2);
		}
		double norm = sqrt(s);
		double tau = (norm*eps+acc)*sqrt(h_new/(b-a));
		/* Såfremt det var succes, så gør jeg klar til næste trin: */
		//printf("Tau is: %g\n",tau);
		//printf("E_i is: %g\n",ei);
		if (ei<tau){
			//printf("Vi går til næste trin\n");
			a += h_new;
			h_new = h_new*pow(tau/ei,0.25)*0.95;
			//printf("Ændrer h: h = %g\n",h_new);
			gsl_vector_memcpy(ya,yh);
			if (n == 2){
				fprintf(outstream,"%g %g %g\n",a,gsl_vector_get(ya,0),gsl_vector_get(ya,1));
			}
			if (n == 3){
				fprintf(outstream,"%g %g %g %g\n",a,gsl_vector_get(ya,0),gsl_vector_get(ya,1),gsl_vector_get(ya,2));
			}
		}
		else {
			//printf("Fejlen var for stor\n");
			h_new = h_new/2;
			//printf("Ændrer h: h = %g\n",h_new);
		}
		
		/* Tjekker lige om vi er færdige: */
		//printf("b - a er lig: %g\n",b-a);
		if ( (b-a) < 0.001) {
			//printf("Færdig\n");
			done = 1;
			gsl_vector_memcpy(yb,ya); /* Kopierer det endelige resultat*/
		}
	}
	gsl_vector_free(yh);
	gsl_vector_free(err);
}

