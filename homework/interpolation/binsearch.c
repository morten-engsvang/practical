#include<gsl/gsl_vector.h>
#include<assert.h>

int binsearch(int n, gsl_vector* x, double z){/* locates the interval for z by bisection */
        //x er vektoren el. listen jeg vil tjekke
        //z er tallet jeg vil finde
        //n er længden af vektoren x
        assert(n>1 && gsl_vector_get(x,0)<=z && z<=gsl_vector_get(x,n-1));//tjekker input
        int i=0, j=n-1;
        while(j-i>1){ //Så længe forskellen er større end én.
                int mid=(i+j)/2; //Giver altid integer (sjovt nok)
                if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
                }
        return i; //Den giver altid punktet der kommer før z :)
        }

