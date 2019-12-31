#include <math.h>
#include <stdio.h>

double gamdev(int ia, long *idum)
/* only for integer paramater ia */
{
	double ran2(long *idum);
	/*void nrerror(char error_text[]);*/
	int j;
	double am, e, s, v1, v2, x, y;
	
	if (ia < 1) {
		printf("Error in routine gamdev");
		} /*nrerror("Error in routine gamdev");*/
	if (ia < 6) { /* use direct method, adding waiting times */
		x=1.0;
		for (j=1;j<=ia;j++) x *= ran2(idum);
		x = -log(x);
	} else {
		do {
			do {
				do {
					v1=ran2(idum);
					v2=2.0*ran2(idum)-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (ran2(idum) > e);
	}
	return x;
}

