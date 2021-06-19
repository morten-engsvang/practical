#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_blas.h>

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
        //Løser Rx=Q^T*b ved back-sub
        backsub(R,x);
	//Output er x som er løsningen til det system
}


void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
        //Udnytter at Q er ortogonal, dvs. transponering giver den inverse
        //Jeg skal derfor bare finde den inverse af den øvre triangulære matrix R
	//Outputtet er B som er den inverse matrix af A.
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

