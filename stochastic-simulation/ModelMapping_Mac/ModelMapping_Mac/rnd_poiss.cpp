#include <math.h>
#define PI 3.141592654

double poidev( double xm, long *idum )
{
	double gammln( double xx );
	double ran2( long *idum );
	static double sq, alxm, g, oldm=(-1.0);
	double em, t, y;

	if ( xm < 12.0 ) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t = 1.0;
		do {
			++em;
			t *= ran2(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm = xm;
			sq = sqrt(2.0*xm);
			alxm = log(xm);
			g = xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y = tan(PI*ran2(idum));
				em = sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
			/* The factor 0.9 is for convenience */
		} while (ran2(idum) > t);
	}
	return em;
}