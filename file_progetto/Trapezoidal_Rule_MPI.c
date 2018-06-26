#include <mpi.h>             /* MPI header file */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double niter = 1E7, d = 1E-7, d2 = 1E-14; 	/* Parameters for estimation of pi (d e d2), # of iterations default (niter) */

double trapezoidal_rule(int start,int finish){
	int i;
	double result = 0.0, x2 = 0.0;
	for(i=start; i< finish; i++){
		x2 = d2 * i * i;
		result += 1.0 / (1.0 + x2);
	}
	return result;
}


int main (int argc, char *argv[])
{
	int my_rank;  		/* Rank current process */
  	int rank;    		/* Loop variable for the processes */
  	int num_proc;   	/* Total number of processes */
  	
  	int quotient;  		/* # of iterations minimum for each process: niter / num_proc */
  	int rem;  		/* Define if someone process take an extra iteration: niter % num_proc */
  	int sub_iter;    	/* Total # of iterations for current process */
  	int sub_start;  	/* Iteration starts for current process */
  

	double pi = 0.0;	/* Estimate of the value of PI */
	double result = 0.0, local_result = 0.0;	/* Global and local result */ 
	double begin, end; 	/* Time elapsed*/

  	/************ SETUP MPI ************/
  	MPI_Init(&argc, &argv);     
  	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);		/* Determine rank of current process */
  	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);	/* Determine # of processes */
  	MPI_Status   status;  				/* Status for receive */ 

	/*********** USAGE ***********/
	if(argc > 1){ 
		if(sscanf(argv[1], "%lf", &niter) != 1){
			if(my_rank == 0)fprintf(stderr, "Invalid Format\nUsage: %s number_of_iterations [Optional - Default 1E7]\n", argv[0]); 
			MPI_Finalize();              	/* Close up MPI  */
			return -1;
		}
		d = pow(10, (int)log10(niter)*-1);
		d2 = pow(10, (int)log10(niter)*-2);
	}
	
	/************ START ************/
  	begin = MPI_Wtime();
	quotient = (int) niter / num_proc;
    	rem = (int) niter % num_proc;
	
	/* For each process calculate where iteration starts and the # of iterations */
	sub_iter = my_rank< rem ? quotient+1 : quotient;
    	if(rem == 0) sub_start = my_rank*quotient;
    	else sub_start = my_rank*quotient+my_rank-(my_rank>rem);

	local_result = trapezoidal_rule(sub_start,sub_start+sub_iter); 	/* current process calculate is local_result */	
    	
	/* Current process send is local_result to Master process (rank = 0) */
	if(my_rank != 0) MPI_Send(&local_result,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	else{ // Master process
		result = local_result;

		/* Get back the local_result from the others */
        	for (rank=1;rank<num_proc;++rank){
             		MPI_Recv(&local_result,1,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD, &status);
             		result += local_result;
        	}
		/* Calculation of the approximation of PI */
        	pi = 4 * d * result;
        	end = MPI_Wtime();
		/************ END ************/
		printf("# iterations = %lf\n", niter);
		printf("# processors = %d\n", num_proc);
		printf("Time execution = %lf\n", end-begin);
		printf("Estimate PI = %.16lf\n", pi);
		printf("Real PI = %.16lf\n", M_PI);
		printf("Error = %.16lf\n", fabs(M_PI - pi));
	}

	MPI_Finalize();              /* Close up MPI  */
	return 0;
}
