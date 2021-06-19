#include<math.h>
#include<complex.h>
#include<stdio.h>

int main(void){
	double G=tgamma(5);
	double B=jn(1,0.5);
	printf("Gamma(5) = %g and Bessel(0.5) = %g\n", G, B);
	complex a=csqrt(-2);
	complex b=cexp(I*M_PI);
	complex c=cexp(I);
	complex d=cpow(I,M_E);
	complex e=cpow(I,I);
	printf("sqrt(-2) = %g + I%g, e^i*pi = %g + I%g , e^i = %g + I%g, i^e = %g + I%g, i^i = %g + I%g\n", creal(a),cimag(a),creal(b),cimag(b),creal(c),cimag(c),creal(d),cimag(d),creal(e),cimag(e));
	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	printf("Float = %.25g, double = %.25lg, long double = %.25Lg\n",x_float,x_double,x_long_double);

return 0;
}
