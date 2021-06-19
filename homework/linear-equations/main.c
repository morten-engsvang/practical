#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_blas.h>
#define RND (double)rand()/RAND_MAX

void backsub(gsl_matrix* U, gsl_vector* c){
	for(int i=c->size-1; i>=0;i--){
		double s = gsl_vector_get(c,i);
		for(int k=i+1;k<c->size;k++){
			s-=gsl_matrix_get(U,i,k)*gsl_vector_get(c,k);
		}
		gsl_vector_set(c,i,s/gsl_matrix_get(U,i,i));
	}
}


void GS_decomp(gsl_matrix* A, gsl_matrix* R){
	//Foretager GR-decomp via. modificeret Gram-Schmidt ortog.
	//Input: n*m (n>=m) matrix A, R er en m*m matrix (tom?)
	//Output: Beregnet R matrix, A er erstatter med Q matricen.
	gsl_vector* ai = gsl_vector_alloc(A->size1);
	gsl_vector* aj = gsl_vector_alloc(A->size1);
	gsl_vector* temp = gsl_vector_alloc(A->size1);	
	//Midlertidig vektor til mellemregninger af samme længde som søjlerne
	for(int i=0;i<A->size2;i++){
		//Kopierer søjle i over i temp1.
		gsl_matrix_get_col(ai,A,i); //Kopierer søjle i over i ai
		//Bruger temp1 til at finde længden af ai
		double square_ai;
		gsl_blas_ddot(ai,ai,&square_ai);
		double ai_length = sqrt(square_ai);
		//Diagonalelementerne (i,i) af R er lig længden af søjlevektorerne i A:
		gsl_matrix_set(R,i,i,ai_length);
		//Normerer min udtagne søjlevektor:
		gsl_vector_scale(ai,1/ai_length);
		//Opdaterer A-matricen
		gsl_matrix_set_col(A,i,ai);
		for(int j=i+1;j<A->size2;j++){
			//Kopierer søjle j over i aj.
			gsl_matrix_get_col(aj,A,j);
			//De resterende (i,j) elementer af R er givet ved følgende, hvor ai nu
			//er normeret pga. tidligere operation
			double R_ij;
			gsl_blas_ddot(ai,aj,&R_ij);
			//Den er dog øvre triangulær:
			gsl_matrix_set(R,i,j,R_ij);
			gsl_matrix_set(R,j,i,0);
			//Jeg opdaterer nu søjlevektor aj
			gsl_vector_memcpy(temp,ai);
			//double square_aj;
			//gsl_blas_ddot(aj,aj,&square_aj);
			//double aj_length = sqrt(square_aj);
			gsl_vector_scale(temp,-1*R_ij);
			gsl_vector_add(aj,temp);
			gsl_matrix_set_col(A,j,aj);
			}
	}
	gsl_vector_free(ai);
	gsl_vector_free(aj);
	gsl_vector_free(temp);
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x); //Danner x ved Q^T*b
	//printf("Debug\n");
	//gsl_vector_fprintf(stdout,x,"%g");
	//Løser Rx=Q^T*b ved back-sub
	backsub(R,x);
	//printf("Debug\n");
	//gsl_vector_fprintf(stdout,x,"%g");
}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
        //Udnytter at Q er ortogonal, dvs. transponering giver den inverse
        //Jeg skal derfor bare finde den inverse af den øvre triangulære matrix R
	gsl_vector* x = gsl_vector_alloc(B->size1);
        gsl_vector* b = gsl_vector_alloc(B->size1);
        for (int i = 0; i<B->size1;i++){
                gsl_vector_set_basis(b,i); //Giver mig i'te basisvektor
                GS_solve(Q,R,b,x);
		gsl_matrix_set_col(B,i,x);
        }
        gsl_vector_free(x);
        gsl_vector_free(b);
}


int main() {
	int n = 3,m = 5;
	printf("Jeg bruger gsl_matrix_fprintf, som printer hvert element linje for linje.\n");
	printf("A-Del1------------------------------------------------------\n");
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_matrix* Q = gsl_matrix_alloc(n,m);
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	gsl_matrix* QTQ = gsl_matrix_alloc(m,m);
	gsl_matrix* QR = gsl_matrix_alloc(n,m);
	//printf("Matricer er allokeret\n");
	for(int i=0;i<A->size1;i++){
		//printf("i=%d\n",i);
		for(int j=0;j<A->size2;j++){
			//printf("j=%d\n",j);
			gsl_matrix_set(A,i,j,RND);
		}
	}
	//printf("Tjek før\n");
	printf("QR-faktorisering via Gram-Schmidt orthogonalisation\n");
	printf("Tilfældig matrix (%d,%d)A =\n",n,m);
	gsl_matrix_fprintf(stdout,A,"%g");
	gsl_matrix_memcpy(Q,A);
	GS_decomp(Q,R);
	//printf("Tjek efter\n");
	printf("Q=\n");
	gsl_matrix_fprintf(stdout,Q,"%g");
	printf("R=(skulle gerne være øvre triangulær)\n");
	gsl_matrix_fprintf(stdout,R,"%g");
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Q,R,0,QR);
	printf("Q*R= (skulle være identisk med A)\n");
	gsl_matrix_fprintf(stdout,QR,"%g");
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Q,Q,0,QTQ);
	printf("Q^T*Q= (skulle give identitetsmatricen)\n");
	gsl_matrix_fprintf(stdout,QTQ,"%g");
	gsl_matrix_free(A);
        gsl_matrix_free(Q);
        gsl_matrix_free(R);
        gsl_matrix_free(QTQ);
        gsl_matrix_free(QR);
	
	printf("A-Del2------------------------------------------------------\n");
	n = 2;
	gsl_matrix* A_square = gsl_matrix_alloc(n,n);
        gsl_matrix* Q_square = gsl_matrix_alloc(n,n);
        gsl_matrix* R_square = gsl_matrix_alloc(n,n);
	gsl_vector* b = gsl_vector_alloc(n);
        gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* tjek = gsl_vector_alloc(n);
	for(int i=0;i<A_square->size1;i++){
                //printf("i=%d\n",i);
                for(int j=0;j<A_square->size2;j++){
                        //printf("j=%d\n",j);
                        gsl_matrix_set(A_square,i,j,RND);
                }
                gsl_vector_set(b,i,RND);
        }
	gsl_matrix_memcpy(Q_square,A_square);
	GS_decomp(Q_square,R_square);
	printf("Tjekker nu om jeg faktisk kan løse et system:\n");
	printf("Tilfældig matrix (%d,%d) A =\n",n,n);
        gsl_matrix_fprintf(stdout,A_square,"%g");
	printf("b=\n");
	gsl_vector_fprintf(stdout,b,"%g");
	GS_solve(Q_square,R_square,b,x);
	printf("Mit fundne x giver Ax=\n");
	gsl_blas_dgemv(CblasNoTrans,1,A_square,x,0,tjek);
	gsl_vector_fprintf(stdout,tjek,"%g");

	printf("B-Del1------------------------------------------------------\n");
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	printf("Jeg finder den inverse af den kvadratiske A matrix fra A-Del2\n");	
	gsl_matrix* inv_tjek = gsl_matrix_alloc(n,n);
	GS_inverse(Q_square,R_square,B);
	printf("Dvs. jeg ser på A=\n");
	gsl_matrix_fprintf(stdout,A_square,"%g");
	printf("Den inverse af A er B=\n");
	gsl_matrix_fprintf(stdout,B,"%g");
	printf("A*B=\n");
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A_square,B,0,inv_tjek);
	gsl_matrix_fprintf(stdout,inv_tjek,"%g");
	printf("B*A=\n");
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,A_square,0,inv_tjek);
	gsl_matrix_fprintf(stdout,inv_tjek,"%g");
	gsl_matrix_free(A_square);
	gsl_matrix_free(Q_square);
	gsl_matrix_free(R_square);
	gsl_matrix_free(B);
	gsl_matrix_free(inv_tjek);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(tjek);

return 0;
}
