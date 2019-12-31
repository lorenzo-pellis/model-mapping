#include <math.h>
#define PI 3.141592654

double bnldev(double pp, int n, long *idum)
/* Returns a floating-point of an integer value (I don't know why) */
{
	double gammln(double xx);
	double ran2(long *idum);
	int j;
	static int nold=(-1);
	double am, em, g, angle, p, bnl, sq, t, y;
	static double pold=(-1.0), pc, plog, pclog, en, oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp); /* Trick to have always p<=0.5 */
	am=n*p; /* Mean value */;
	if (n < 25) { /*Use direct method if nis small */
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran2(idum) < p) ++bnl;
	} else if (am < 1.0) { /* If less than 1 event out of 25, use the Poisson */
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran2(idum);
			if (t <g ) break;
		}
		bnl=(j <= n ? j : n);
	} else { /* Rejection method */
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran2(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran2(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
