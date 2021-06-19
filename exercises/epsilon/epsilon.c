#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(void){
	/*
	//Udfører først 1i
	printf("INT_MAX = %i\n",INT_MAX);
	//Udfører 1i for while-loop
	int i = 1;
	while (i+1>i) {
		i++;
	}
	printf("My max int via. while loop = %i\n",i);
	//Udfører 1i for for-loop
	for (i=1;i+1>i;i++){
	}
	printf("My max int via. for loop = %i\n",i);

	//Udfører 1i for do while loop
	i = 1;
	do{i++;}
	while(i+1>i);
	printf("My max int via. do-while loop = %i\n",i);

	//Udfører derefter 1ii
	printf("INT_MIN = %i\n",INT_MIN);
	i = -1;
	while (i-1<i) {
		i--;
	}
	printf("My min int via. while loop = %i\n",i);

	for (i=-1;i-1<i;i--) {
	}
	printf("My min int via. for loop = %i\n",i);

	i = -1;
	do{i--;}
	while(i-1<i);
	printf("My min int via. do-while loop = %i\n",i);
	*/

	//Udfører til sidst 1iii
	//Først for float
	
	printf("FLT_EPSILON = %f\n",FLT_EPSILON);
	float x = 1;
	while (1+x!=1) {
		x/=2;
	}
	x *= 2;
	printf("My float epsilon via. while loop = %f\n",x);
	for (x=1;1+x!=1;x/=2) {
	}
	x *= 2;
	printf("My float epsilon via. for loop = %f\n",x);
	x = 1;
	do {x/=2;}
	while (1+x!=1);
	x *= 2;
	printf("My float epsilon via. do-while loop = %f\n",x);
	
	//Så for double
	printf("DBL_EPSILON = %g\n",DBL_EPSILON);
	double x1 = 1;
	while (1+x1!=1) {
		x1/=2;
	}	
	x1 *= 2; //Det ønskede epsilon er det sidste hvor betingelsen er true, vi går ét skridt for langt.
	printf("My double epsilon via. while loop = %g\n",x1);
	for (x1=1;1+x1!=1;x1/=2){
	}
	x1 *= 2;
	printf("My double epsilon via. for loop = %g\n",x1);
	
	x1 = 1;
	do {x1/=2;}
	while (1+x1!=1);
	x1 *= 2;
	printf("My double epsilon via. do-while loop = %g\n",x1);
	
	printf("LDBL_EPSILON = %Lg\n",LDBL_EPSILON);
	long double x2 = 1;
	while (1+x2!=1) {
		x2/=2;
	}
	x2 *= 2;
	printf("My long double epsilon via. while loop = %Lg\n",x2);
	
	for (x2=1;1+x2!=1;x2/=2) {
	}
	x2 *= 2;
	printf("My long double epsilon via. for loop = %Lg\n",x2);
	
	x2 = 1;
	do {x2/=2;}
	while (1+x2!=1);
	x2 *= 2;
	printf("My long double epsilon via. do-while loop = %Lg\n",x2);
	
	//Udfører nu for 2i
	int max=INT_MAX/2;
	float sum_up_float = 0;
	for (int i=1;i<=max;i++) {
		sum_up_float += 1.0f/i;
	}
	float sum_down_float = 0;
	for (int i=max;i>=1;i--) {
		sum_down_float += 1.0f/i;
	}
	printf("sum_up_float = %f\n",sum_up_float);
	printf("sum_down_float = %f\n",sum_down_float);
	printf("Forskellen skyldes at vi bruger float typen. Den har en begrænset præcision\nJeg forsøger derfor med max = int_max/1.5\nDet burde dog ikke virke da det giver led der er næsten lig 0\n");
        max = INT_MAX/1.5f;
	sum_up_float = 0;
	for (int i=1;i<=max;i++) {
		sum_up_float += 1.0f/i;
	}
	sum_down_float = 0;
	for (int i=max;i>=1;i--) {
		sum_down_float += 1.0f/i;
	}
	printf("sum_up_float = %f\n",sum_up_float);
	printf("sum_down_float = %f\n",sum_down_float);
	printf("Jeg forsøger nu det samme med double\n");
	max=INT_MAX/2;
	double sum_up_double = 0;
	for (int i=1;i<=max;i++) {
		sum_up_double += 1.0/i;
	}
	double sum_down_double = 0;
	for (int i=max;i>=1;i--) {
		sum_down_double += 1.0/i;
	}
	printf("sum_up_double = %g\n",sum_up_double);
	printf("sum_down_double = %g\n",sum_down_double);
	printf("De to summer konvergerer, da vi nu har mere præcision på vores tal pga. datatype\n");
return 0;
}
