#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<math.h>

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
