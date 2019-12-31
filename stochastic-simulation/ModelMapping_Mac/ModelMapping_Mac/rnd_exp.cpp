#include <math.h>
#include <stdio.h>

double expdev(long *idum)
/* this is exp with mean 1. For b*exp(-b*x) the output must be divided by b */
{
	double ran2(long *idum);
	double dum;

	do
		dum=ran2(idum);
	while (dum == 0.0);
	return -log(dum);
}
