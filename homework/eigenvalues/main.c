#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<math.h>
#define RND (double)rand()/RAND_MAX

//Vi skal i sidste endnu udføre en klassisk Jacobi transformation
//hvor række efter række elimineres.

void vector_print(const char* s, gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%9.3f",gsl_vector_get(v,i));
	printf("\n");
}

void matrix_print(const char* s, gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++)printf("%9.3f",gsl_matrix_get(A,i,j));
		printf("\n");
	}
}



void timesJ(gsl_matrix* A, int p, int q, double theta){
	//Input: Matrix A, der ønskes ganget med Jacobi matricen defineret ved J(p,q,theta)
	//Jacobi matricen skal fjerne element A_pq, definition givet i eigenvalue kapitlet.
	//Output: Matrix A, der nu er ganget med J fra højre
	//Starter ud med at definere forkortelser for cos(theta) og sin(theta)
	double c = cos(theta), s = sin(theta);
	//Jeg finder nu elementerne A_ip og A_iq efter multiplikation
	//Hvor i betegner rækkerne, da jeg ender med at gange ind på alle rækkerne
	for(int i=0;i<A->size1;i++){
		//Først multiplikation med den p'te søjle i Jacobi matricen
		//som opdaterer den p'te søjle i A ud fra p'te og q'te element i række i.
		double A_ip = c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		//Derefter multiplikation med den q'te søjle i Jacobi matricen
		double A_iq = s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		//Jeg opdaterer nu min række i A matricen:
		gsl_matrix_set(A,i,p,A_ip);
		gsl_matrix_set(A,i,q,A_iq);
	}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	//Input: Matrix A, der ønskes ganget med Jacobi matricen defineret ved J(p,q,theta)
	//Output: Matrix A, der nu er ganget med J fra venstre
	double c = cos(theta), s = sin(theta);
	//Multiplikation fra den anden side medfører at jeg skal bruge j til at definere søjler.
	//og vi ender med at gange ind på alle søjlerne i A.
	for(int j=0;j<A->size2;j++){
		//Først multiplikation med den p'te række i Jacobi matricen
		//som opdaterer den p'te række i A ud fra p'te og q'te element i søjle j.
		double A_pj = c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		//Så med den q'te række i Jacobi matricen:
		double A_qj = -s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		//Til sidst opdateres søjlen i A matricen:
		gsl_matrix_set(A,p,j,A_pj);
		gsl_matrix_set(A,q,j,A_qj);
	}
}

void jacobi_diag(gsl_matrix* A, gsl_matrix* V){
	//Input: Matrix A der skal undersøges, matrix V der er enhedsmatricen.
	//Jeg laver loops, der går gennem alle de mulige (p,q) og udfører Jacobi transformationen
	//Output: Diagonaliseret A matrix, V matrix af egenvektorer.
	//A dannes ved J^T*A*J, og V dannes ved I*J
	int change = 1;
	double delta = 0.000001;
	while(change!=0){
		change = 0;
		for(int p=0;p<A->size1;p++){
			for(int q=p+1;q<A->size1;q++){
				//Trækker først elementer ud til at tjekke om jeg skal udføre
				//transformationen (og danne den theta jeg skal bruge)
				double A_pq = gsl_matrix_get(A,p,q);
				double A_pp = gsl_matrix_get(A,p,p);
				double A_qq = gsl_matrix_get(A,q,q);
				//Der er fejl i ligning elleve, vi skal også gange med en halv
				//for at få theta
				double theta = 0.5*atan2(2*A_pq,A_qq-A_pp);
				double c = cos(theta), s = sin(theta);
				//Jeg finder udvalgte elementer efter transformation
				//Når vi er færdige vil de være egenværdierne:
				double A_pp_updated = c*c*A_pp-2*s*c*A_pq+s*s*A_qq;
				double A_qq_updated = s*s*A_pp+2*s*c*A_pq+c*c*A_qq;
				//Jeg kan så tjekke om der sker en ændring ved transformation
				//for at se om jeg skal ændre flag til at fortsætte:
				//printf("Debug: forskellene er: %g og %g\n",fabs(A_pp_updated-A_pp), fabs(A_qq_updated-A_qq));
				if((fabs(A_pp_updated-A_pp))>delta || (fabs(A_qq_updated-A_qq))>delta){
					//Opdater flag, da vi ser en ændring:
					change = 1;
					timesJ(A,p,q,theta); //A*J
					Jtimes(A,p,q,-theta); //J^T*A
					//Husk at opdatere vektorer:
					timesJ(V,p,q,theta);
				}
			}
		}
	}
}


int main(void){
	//Skal lave en tilfældig reel symmetrisk matrix A
	//Udføre eigenværdi decomp
	int n = 3;
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* D = gsl_matrix_alloc(n,n);
	gsl_matrix* I = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix* tjek1 = gsl_matrix_alloc(n,n);
	gsl_matrix* tjek2 = gsl_matrix_alloc(n,n);
	gsl_vector* b = gsl_vector_alloc(n);
	for(int i = 0;i<n;i++){
		gsl_vector_set_basis(b,i);
		gsl_matrix_set_col(I,i,b);
		for(int j = 0;j<n;j++){
			double number = RND;
			gsl_matrix_set(A,i,j,number);
			gsl_matrix_set(A,j,i,number);
		}
	}
	//Jeg kan nu udføre min jacobi transformation:
	gsl_matrix_memcpy(D,A);
	gsl_matrix_memcpy(V,I);
	jacobi_diag(D,V);
	printf("DelA---------------------------------------------\n");
	matrix_print("V=",V);
	printf("Tester V^T*A*V=D, printer de to sider:\n");
	matrix_print("D=",D);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,A,0,tjek1);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,tjek1,V,0,tjek2);
	matrix_print("V^T*A*V=",tjek2);
	printf("Tester nu V*D*V^T==A, printer de to sider:\n ");
	matrix_print("A=",A);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,D,0,tjek1);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,tjek1,V,0,tjek2);
	matrix_print("V*D*V^T=",tjek2);
	printf("Til sidst tjekker jeg om V^T*T=I:");
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,V,0,tjek2);
	matrix_print("V^T*T=",tjek2);

	gsl_matrix_free(A);
	gsl_matrix_free(D);
	gsl_matrix_free(I);
	gsl_matrix_free(V);
	gsl_matrix_free(tjek1);
	gsl_matrix_free(tjek2);
	gsl_vector_free(b);
	printf("DelB---------------------------------------------\n");
	//Bygger Hamilton operatoren:
	printf("Jeg har konstrueret og diagonaliseret Hamilton operatoren.\n");
	printf("Egenværdierne kan ses sammen med de analytiske værdier i energy.txt\n");
	printf("Egenfunktionerne kan ses i eigen.png og egenværdierne er plottet i energy.png");
	n=50;
	double s=1.0/(n+1);
	gsl_matrix* H = gsl_matrix_alloc(n,n);
	for(int i=0;i<n-1;i++){
		gsl_matrix_set(H,i,i,-2);
		gsl_matrix_set(H,i,i+1,1);
		gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,n-1,n-1,-2);
	gsl_matrix_scale(H,-1/s/s);
	//Diagonalisér matrixen:
	gsl_matrix* X = gsl_matrix_alloc(n,n);
	gsl_matrix_set_identity(X);
	jacobi_diag(H,X);
	//Tjekker om det er korrekt:
	FILE * energy = fopen("energy.txt","w");
	//fprintf(energy,"k exact calculated\n");
	for (int k=0; k < n/3; k++){
		double exact = M_PI*M_PI*(k+1)*(k+1);
		double calculated = gsl_matrix_get(H,k,k);
		fprintf(energy,"%i %g %g\n",k,calculated,exact);
	}
	fclose(energy);
	FILE * eigen = fopen("eigen.txt","w");
	for(int k=0;k<3;k++){
		fprintf(eigen,"#INDEX %i\n",k);
		fprintf(eigen,"0 0 0\n");
		for(int i=0;i<n;i++){
			double exact = sqrt(2)*sin((k+1)*M_PI*(i+1.0)/(n+1));
			double calc = gsl_matrix_get(X,i,k)*(1/sqrt(s))*pow(-1,k); //Den normaliserede funktion. Hver anden har fået omvendt fortegn, er ikke sikker på hvorfor, men må gerne gange med en konstant.
			fprintf(eigen,"%g %g %g\n",(i+1.0)/(n+1), calc, exact);
		}
		fprintf(eigen,"1 0 0\n\n\n");
	}	
	fclose(eigen);
	
	
return 0;
}
