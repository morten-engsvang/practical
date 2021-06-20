#include<gsl/gsl_vector.h>
#include<math.h>

double open4point(double f(double), double a, double b, double delta, double eps, double f2, double f3){
	/* Beregner rekursivt for 4 punkter i et åbent interval fra a til b */
	double f1 = f(a+(1.0/6)*(b-a)), f4 = f(a+(5.0/6)*(b-a));
	/* Jeg kan nu beregne mine approximative integraler */
	double Q = (2.0*f1+f2+f3+2.0*f4)*(b-a)/6; //Trapez reglen, sum af vægt ganget med funktionsværdi
	double q = (f1+f2+f3+f4)*(b-a)/6; //Rektangel reglen
	double tolerance = delta+eps*fabs(Q); //Den tilladte fejl
	double error = fabs(Q-q);
	/*
	if(limit==0){
		fprintf(stderr,"Recursion limit reached\n");
		return Q;
	}
	*/
	
	if(error < tolerance){
		return Q;
	}
	else {
		double Q1 = open4point(f,a,(a+b)/2,delta/sqrt(2),eps,f1,f2);
		double Q2 = open4point(f,(a+b)/2,b,delta/sqrt(2),eps,f1,f2);
		return Q1+Q2;
	}
}

double recursive_integrate(double f(double), double a, double b, double delta, double eps){
	/* Vi får givet ende punkterne og kan beregne funktionsværdierne ved punkt 3 og 4 fra ligning 48, pga. måden open4point er sat op */
	double f2 = f(a+(2.0/6)*(b-a)), f3 = f(a+(4.0/6)*(b-a));
	/* Laver en rekursionsgrænse */
	/* Starter rekursionen */
	return open4point(f,a,b,delta,eps,f2,f3);
}

double clenshaw_curtis(double f(double),double a,double b,double delta,double eps){
	/* Clenshaw-Curtis variabel transformation af et generelt integral fra a til b */
	double g(double y){
		return f( (a+b)/2+(b-a)/2*cos(y) )*sin(y)*(b-a)/2;
	}
	return recursive_integrate(g,0,M_PI,delta,eps);
}
