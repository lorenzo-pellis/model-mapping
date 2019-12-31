#define r2_IM1 2147483563
#define r2_IM2 2147483399
#define r2_AM (1.0/r2_IM1)
#define r2_IMM1 (r2_IM1-1)
#define r2_IA1 40014
#define r2_IA2 40692
#define r2_IQ1 53668
#define r2_IQ2 52774
#define r2_IR1 12211
#define r2_IR2 3791
#define r2_NTAB 32
#define r2_NDIV (1+r2_IMM1/r2_NTAB)
#define r2_EPS 1.2e-7
#define r2_RNMX (1.0-r2_EPS)

double ran2(long *idum)
/* idum must be a negative integer */
/* do not alter it after the first time */
{
	int j;
	long k;
	static long r2_idum2=123456789;
	static long r2_iy=0;
	static long r2_iv[r2_NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		r2_idum2=(*idum);
		for (j=r2_NTAB+7;j>=0;j--) {
			k=(*idum)/r2_IQ1;
			*idum=r2_IA1*(*idum-k*r2_IQ1)-k*r2_IR1;
			if (*idum < 0) *idum += r2_IM1;
			if (j < r2_NTAB) r2_iv[j] = *idum;
		}
		r2_iy=r2_iv[0];
	}
	k=(*idum)/r2_IQ1;
	*idum=r2_IA1*(*idum-k*r2_IQ1)-k*r2_IR1;
	if (*idum < 0) *idum += r2_IM1;
	k=r2_idum2/r2_IQ2;
	r2_idum2=r2_IA2*(r2_idum2-k*r2_IQ2)-k*r2_IR2;
	if (r2_idum2 < 0) r2_idum2 += r2_IM2;
	j=r2_iy/r2_NDIV;
	r2_iy=r2_iv[j]-r2_idum2;
	r2_iv[j] = *idum;
	if (r2_iy < 1) r2_iy += r2_IMM1;
	if ((temp=r2_AM*r2_iy) > r2_RNMX) return r2_RNMX;
	else return temp;
}