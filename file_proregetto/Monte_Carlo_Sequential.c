#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define SEED 35791246

int main (int argc, char *argv[])
{
	int niter = 0;
	double x, y;
	int i, count = 0; /* # of points in the 1st quadrant of unit circle */
	double z;
	double pi;
	printf("Enter the number of iterations used to estimate pi: ");
	scanf("%d", &niter);

	clock_t begin = clock(), end;

	/* nitialize randome numbers */
	srand(SEED);
	count = 0;
	
	for(i=0; i<niter; i++){
		x = (double) rand() / RAND_MAX;
		y = (double) rand() / RAND_MAX;	
		z = x*x + y*y;
		if(z<=1) count++;	
	}
	pi = (double) count / niter * 4;
	end = clock(); 
	printf(" # of trials= %d, estimate of pi is %g\nTime=%lf\n", niter, pi,((double)(end-begin))/CLOCKS_PER_SEC);
	return 0;
}
