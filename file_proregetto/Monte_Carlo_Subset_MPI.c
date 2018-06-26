#include <mpi.h> 	/* MPI header file */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define SEED 35791246

int main (int argc, char *argv[])
{
	double x, y, z, pi;			/* Pseudorandom point (x,y), test point (z), estimate pi*/
	int i, count = 0, local_count = 0; 	/* # of points in the 1st quadrant of unit circle */
	int niter;				/* # of iterations */
	
	int my_rank;  				/* Rank Process */
  	int rank;    				/* Loop variable for the processes */
  	int num_proc;   			/* Total number of processes */

	double begin, end; 			/* Time elapsed*/

	int quotient;  				/* # of iterations minimum for each process: N / num_proc */
  	int rem;  				/* Define if someone process take an extra iteration: N % num_proc */
  	int sub_len;    			/* Total # of iterations for current process */
  			
  	/************ SETUP MPI ************/
  	MPI_Init(&argc, &argv);     
  	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);		/* Determine rank of current process */
  	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);	/* Determine # of processes */
  	MPI_Status   status;  				/* Status for receive */
	
	/*********** USAGE ***********/
	if(argc <= 1){
		if(my_rank==0)fprintf(stderr, "Usage: %s number_of_iterations \n", argv[0]);
		MPI_Finalize();
		return -1;
	}

  	if(sscanf(argv[1], "%d", &niter) != 1) {
		if(my_rank==0) fprintf(stderr, "Invalid Format (1st par must be int) \nUsage: %s number_of_iterations \n", argv[0]); 
		MPI_Finalize();
		return -1;
	}

	/************ START ************/
  	begin = MPI_Wtime();
	srand(SEED);					/* Initialize random numbers */
	
	quotient = (int) niter / num_proc;
    	rem = (int) niter % num_proc;

	/* For each process calculate the # of iterations */
	sub_len = my_rank< rem ? quotient+1 : quotient;

	for(i=0; i<sub_len; i++){
		x = (double) rand() / RAND_MAX;
		y = (double) rand() / RAND_MAX;	
		z = x*x + y*y;
		if(z<=1) local_count++;	
	}

	/* Current process send is count to Master process (rank = 0) */
	if(my_rank != 0) MPI_Send(&local_count,1,MPI_INT,0,0,MPI_COMM_WORLD); 
	else{ // Master process
		count = local_count;

		/* Get back the local_count from the others */
        	for (rank=1;rank<num_proc;++rank){
             		MPI_Recv(&local_count,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD, &status);
             		count += local_count;
        	}
		/* Calculation of the approximation of PI */
		pi = (double) count / niter * 4;
		end = MPI_Wtime();
		/************ END ************/
		printf("# of trials = %d\n", niter);
		printf("# processors = %d\n", num_proc);
		printf("Time execution = %lf\n", end-begin);
		printf("Estimate PI = %.16lf\n", pi);
		printf("Real PI = %.16lf\n", M_PI);
		printf("Error = %.16lf\n", fabs(M_PI - pi));
	}

	MPI_Finalize();              /* Close up MPI */
	return 0;
}
