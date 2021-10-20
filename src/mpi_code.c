#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#define  Max(a,b) ((a)>(b)?(a):(b))

#define  N   ((1 << 6) +2)
double   maxeps = 0.1e-7;
int itmax = 100;
int eps;

int main(int argc, char **argv)
{
    
    MPI_Init(&argc, &argv);

    int rank;
    int size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double time_start = MPI_Wtime();    
    if (size > N - 2) 
	{
        size = N - 2;
    }

    if (rank < size) 
	{
	    int part_size = (N - 2) / size;
	    int ind_first = part_size * rank, ind_last = ind_first + part_size;
    	if (rank == size - 1) 
		{ //last process
        	ind_last = N - 2;
    	} 
    	int *eps_achieved = NULL; // necessary only for zero process 

    	double *matrix = (double *)calloc(N * N * (ind_last - ind_first + 2), sizeof(double));
    	matrix += N * N;
    	int ind1, ind2, ind3;
    	for (ind1 = 0; ind1 < ind_last - ind_first; ind1++) 
		{
        	for (ind2 = 0; ind2 <= N - 1; ind2++) 
			{
            	for (ind3 = 0; ind3 <= N - 1; ind3++) 
				{
                	if (ind2 == 0 || ind3 == 0 || ind3 == N - 1 || ind2 == N - 1 ) {
                    	*(matrix + ind1 * N * N + ind2 * N + ind3) = 0;
   					} else {
   						*(matrix + ind1 * N * N + ind2 * N + ind3) = 4. + (ind_first + ind1 + 1) + ind2 + ind3; 
   					}
   				}
   			}
		}
		int should_calculate = 1;
	
		MPI_Status status;
//matrix initialized
		if (rank != size - 1) 
		{ //for last process next_values is always zero submatrix
        	for (ind2 = 0; ind2 <= N - 1; ind2++) 
			{
            	for (ind3 = 0; ind3 <= N - 1; ind3++) 
				{
                	if (ind2 == 0 || ind3 == 0 || ind3 == N - 1 || ind2 == N - 1 ) {
                    	*(matrix + N * N * (ind_last - ind_first)+ ind2 * N + ind3) = 0;
                	} else {
                    	*(matrix + N * N * (ind_last - ind_first) + ind2 * N + ind3) = 4. + ind2 + ind3 + ind_last - ind_first;
                	}
            	}
        	}
    	}

    	if (!rank) 
		{
        	eps_achieved = (int *)calloc(size, sizeof(int));
    	}
    	int it;
    	for (it = 0; it < itmax; it++) 
		{
        	int flag;
        	if (!rank) 
			{
// check, if have messages about eps achievement from other processes
            	for (int process_id = 1; process_id < size; process_id++) 
				{
                	flag = 0;
                	MPI_Iprobe(process_id, itmax * 2, MPI_COMM_WORLD, &flag, &status);
                	if (flag) {
                    	eps_achieved[process_id] = 1;
                    	MPI_Recv(&flag, 1, MPI_INT, process_id, itmax * 2, MPI_COMM_WORLD, &status); // take this
                    //message from queue
                	} 
           	 	}
            	flag = 1;
            	for (int i = 0; i < size; i++) 
				{
                	if (!eps_achieved[i]) 
					{
                    	flag = 0;
                    	break;
                	}
            	}
            	if (flag) 
				{
                	for (int i = 1; i < size; i++) 
					{ 
                    	MPI_Send(&flag, 1, MPI_INT, i, itmax * 2, MPI_COMM_WORLD);
                	}
            	}
        	} 
			else 
			{
            	MPI_Iprobe(0, itmax * 2, MPI_COMM_WORLD, &flag, &status);
            	if (flag) 
				{
                	MPI_Recv(&flag, 1, MPI_INT, 0, itmax * 2, MPI_COMM_WORLD, &status);
                	should_calculate = 0;
           	 	}
        	}
    // get values
        	if (rank) 
			{ // process zero does not recieves data from other processes here
            //on first iteration left values could be counted by process itself
            	MPI_Recv(matrix - N * N, N * N, MPI_DOUBLE, rank - 1, it, MPI_COMM_WORLD, &status);
        	}
        	if (it && rank != size - 1) 
			{ // on first iteration next_values is already correct,
            //for last process next_values is ALWAYS correct
            	MPI_Recv(matrix + N * N * (ind_last - ind_first), N * N, MPI_DOUBLE, rank + 1, it - 1, MPI_COMM_WORLD, &status);
        	}
        // process values
        	eps = 0;
        	if (should_calculate) 
			{
        		for (ind1 = 0; ind1 < ind_last - ind_first; ind1++) 
				{
            		for (ind2 = 1; ind2 < N - 1; ind2++) 
					{
                		for (ind3 = 1; ind3 < N - 1; ind3++)
						{
                    		double old = *(matrix + ind1 * N * N + ind2 * N + ind3);            
                    		*(matrix + ind1 * N * N + ind2 * N + ind3) = (*(matrix + (ind1 - 1) * N * N + ind2 * N + ind3) + *(matrix + ind1 * N * N + (ind2 - 1) * N + ind3) + *(matrix + ind1 * N * N + (ind2 + 1) * N + ind3) + *(matrix + ind1 * N * N + ind2 * N + ind3 - 1) + *(matrix + ind1 * N * N + ind2 * N + ind3 + 1) + *(matrix + (ind1 + 1) * N * N + ind2 * N + ind3)) / 6.;
							eps = Max(eps, fabs(old - *(matrix + ind1 * N * N + ind2 * N + ind3)));
						}
					}
				}
			}
			if (rank != size - 1)
			{
    	        MPI_Send(matrix + (ind_last - ind_first - 1) * N * N, N * N, MPI_DOUBLE, rank + 1, it, MPI_COMM_WORLD);
        	}
        	if (rank && it != itmax - 1)
			{
            	MPI_Send(matrix, N * N, MPI_DOUBLE, rank - 1, it, MPI_COMM_WORLD);
        	}
        	if (should_calculate)
			{
        		if (eps < maxeps)
				{
            // notify zero process about this
            		if (!rank)
					{
                		eps_achieved[0] = 1;
            		}
				   	else
					{
                		int achieved = 1;
                		MPI_Send(&achieved, 1, MPI_INT, 0, itmax * 2, MPI_COMM_WORLD);
            		}
         		}

        	}
		}
//now we should sum all counted values and sent them to process 0
//process 0 should sum all recieved values and print them
    
	    double S = 0.;
   	 	for (ind1 = 0; ind1 < ind_last - ind_first; ind1++)
		{
        	for (ind2 = 1; ind2 < N - 1; ind2++)
			{
            	for (ind3 = 1; ind3 < N - 1; ind3++)
				{
                	S += *(matrix + ind1 * N * N + ind2 * N + ind3) * (ind3 + 1) * (ind2 + 1) * (ind1 + ind_first + 2) / (N * N * N);
            	}
        	}
    	}
    	if (rank)
		{
        	MPI_Send(&S, 1, MPI_DOUBLE, 0, itmax + 1, MPI_COMM_WORLD);
    	}

    	if (!rank)
		{
        	double S_tmp;
        	for (ind1 = 1; ind1 < size; ind1++)
			{
            	MPI_Recv(&S_tmp, 1, MPI_DOUBLE, ind1, itmax + 1, MPI_COMM_WORLD, &status);
            	S += S_tmp;
        	}
    	}
    	if (!rank)
		{
			fprintf(stderr, "S = %lf\n", S);
    	}
    	free(matrix - N * N);
    }
    MPI_Finalize();
	return 0;
}
