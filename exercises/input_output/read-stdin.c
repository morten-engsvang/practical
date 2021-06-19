#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
	double x; 
	int input;
	do {
		input=fscanf(stdin,"%lg",&x);
		printf("x=%g sin(x)=%g cos(x)=%g\n",x,sin(x),cos(x));
	}
	while(input!=EOF);
return 0;
}
