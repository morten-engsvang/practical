#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

//Stjæler vector_print skamløst
void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for (int i=0; i < v->size;i++) {
		printf("%10g ",gsl_vector_get(v,i));
	}
	printf("\n");
}



int main(){
	//Jeg definerer A-matrix, kopi, b-vektor og de tomme x og y vektorer
	gsl_matrix* A = gsl_matrix_alloc(3,3);
	gsl_matrix* Acopy = gsl_matrix_alloc(3,3); //householder ødelægger A
	gsl_vector* b = gsl_vector_alloc(3);
	gsl_vector* x = gsl_vector_alloc(3);
	gsl_vector* y = gsl_vector_alloc(3);
	gsl_matrix_set(A,0,0,6.13);//Den forventer at værdien er en double
	gsl_matrix_set(A,1,0,8.08);
	gsl_matrix_set(A,2,0,-4.36);
	gsl_matrix_set(A,0,1,-2.90);
	gsl_matrix_set(A,1,1,-6.31);
	gsl_matrix_set(A,2,1,1.00);
	gsl_matrix_set(A,0,2,5.86);
	gsl_matrix_set(A,1,2,-3.89);
	gsl_matrix_set(A,2,2,0.19);
	gsl_matrix_memcpy(Acopy,A);
	gsl_vector_set(b,0,6.23);
	gsl_vector_set(b,1,5.37);
	gsl_vector_set(b,2,2.29);
	//Jeg løser Ax=b ved householder solver
	gsl_linalg_HH_solve(A,b,x);
	vector_print("Solutions is x:",x);
	//Tjekker mit resultat
	gsl_blas_dgemv(CblasNoTrans,1,Acopy,x,0,y);
	vector_print("A*x then gives b:",y);
	vector_print("It should be equal to b:",b);
	gsl_matrix_free(A);
	gsl_matrix_free(Acopy);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(y);
return 0;
}
