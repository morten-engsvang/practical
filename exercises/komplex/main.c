#include"komplex.h"
#include<stdio.h>
#define TINY 1e-6

int main(void){
	komplex a = {1,2}, b = {3,4};

	printf("testing komplex add:\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should be = ", R);
	komplex_print("a+b is actually = ", r);
}
