#include<stdio.h>
#include<math.h>

int main (int arc, char *argv[]) {
	FILE* input=fopen(argv[1],"r");
	FILE* output=fopen("out.txt","w");
	int items;
	double x;

	do{
		items=fscanf(input,"%lg",&x);
		fprintf(output, "x=%g sin(x)=%g cos(x)=%g\n",x,sin(x),cos(x));
	}
	while(items!=EOF);


	fclose(input);
	fclose(output);
return 0;
}
