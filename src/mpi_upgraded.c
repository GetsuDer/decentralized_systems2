#include <fcntl.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"
#define  Max(a,b) ((a)>(b)?(a):(b))

#define  N   ((1 << 6) +2)
#define RESTORE_TAG 1
#define EPS_TAG 2
#define END_TAG 3
double   maxeps = 0.1e-7;
int itmax = 100;
int eps;

struct env
{
	int ind_first; // first index to work with
	int ind_last; // last index to work with
	int prev_proc; // rank of process which gives previous values
	int next_proc; // rank of process to which must give values
	int rank; // rank of current process
	int size; // number of existing processes
	int it; // last successfull iteration
	int *live_procs; // array of currently live processes
};

void
fixing_all(struct env *env)
{
	printf("In fixing_all proc %d/%d on it %d\n", env->rank, env->size, env->it);
	while (1)
	{

	}
	return;
}

int main(int argc, char **argv)
{
	srandom(time(NULL));
    int bad_rank = random();
	int bad_it = random() % itmax + 1;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	struct env e;
	int rc;
	int err_shift = 1;
	int rearrange_buff[4];
	MPI_Request request;

    MPI_Comm_size(MPI_COMM_WORLD, &e.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &e.rank);

	if (!(bad_rank % e.size))
	{
		bad_rank++;
	}

	if (!e.rank) 
		printf("bad_it %d bad_rank %... %d\n", bad_it, bad_rank % e.size); 
    if (e.size > N - 2) 
	{
        e.size = N - 2;
    }

    if (e.rank < e.size) 
	{
		e.live_procs = (int *)calloc(e.size, sizeof(int));
		for (int i = 0; i < e.size; i++)
		{
			e.live_procs[i] = 1;
		}
	    int part_size = (N - 2) / e.size;
	    e.ind_first = part_size * e.rank;
	    e.ind_last = e.ind_first + part_size;
		e.prev_proc = e.rank? e.rank - 1 : e.size - 1;
		e.next_proc = (e.rank == e.size - 1) ? 0 : e.rank + 1;
    	if (e.rank == e.size - 1) 
		{ //last process
        	e.ind_last = N - 2;
    	} 
    	int *eps_achieved = NULL; // necessary only for zero process 

    	double *matrix = (double *)calloc(N * N * (e.ind_last - e.ind_first + 2), sizeof(double));
    	matrix += N * N;
    	int ind1, ind2, ind3;
    	for (ind1 = 0; ind1 < e.ind_last - e.ind_first; ind1++) 
		{
        	for (ind2 = 0; ind2 <= N - 1; ind2++) 
			{
            	for (ind3 = 0; ind3 <= N - 1; ind3++) 
				{
                	if (ind2 == 0 || ind3 == 0 || ind3 == N - 1 || ind2 == N - 1 ) {
                    	*(matrix + ind1 * N * N + ind2 * N + ind3) = 0;
   					} else {
   						*(matrix + ind1 * N * N + ind2 * N + ind3) = 4. + (e.ind_first + ind1 + 1) + ind2 + ind3; 
   					}
   				}
   			}
		}
		int should_calculate = 1;
	
		MPI_Status status;
//matrix initialized
		if (e.rank != e.size - 1) 
		{ //for last process next_values is always zero submatrix
        	for (ind2 = 0; ind2 <= N - 1; ind2++) 
			{
            	for (ind3 = 0; ind3 <= N - 1; ind3++) 
				{
                	if (ind2 == 0 || ind3 == 0 || ind3 == N - 1 || ind2 == N - 1 ) {
                    	*(matrix + N * N * (e.ind_last - e.ind_first)+ ind2 * N + ind3) = 0;
                	} else {
                    	*(matrix + N * N * (e.ind_last - e.ind_first) + ind2 * N + ind3) = 4. + ind2 + ind3 + e.ind_last - e.ind_first;
                	}
            	}
        	}
    	}

    	if (!e.rank) 
		{
        	eps_achieved = (int *)calloc(e.size, sizeof(int));
    	}
		e.it = 0;
    	for (; e.it < itmax;) 
		{
			if ((e.it == bad_it) && e.rank == (bad_rank % e.size))
			{
				// raise error
				printf("On iteration %d process %d/%d decided to die.\n", e.it, e.rank, e.size);
				exit(1);
			}

        	int flag;
        	if (!e.rank) 
			{
// check, if have messages about eps achievement from other processes
            	for (int process_id = 1; process_id < e.size; process_id++) 
				{
					if (!e.live_procs[process_id])
					{
						// this process is already dead, forget about it
						continue;
					}
                	flag = 0;
                	rc = MPI_Iprobe(process_id, EPS_TAG + itmax * err_shift, MPI_COMM_WORLD, &flag, &status);
					if (rc != MPI_SUCCESS)
					{
						printf("Process %d on it %d recieved not OK 1\n", e.rank, e.it);
						fixing_all(&e);
						continue; // go to next iteration
					}

                	if (flag) {
                    	eps_achieved[process_id] = 1;
                    	rc = MPI_Recv(&flag, 1, MPI_INT, process_id, EPS_TAG + itmax * err_shift, 
								MPI_COMM_WORLD, &status); 
						// take this
                    	//message from queue
						if (rc != MPI_SUCCESS)
						{
							printf("Process %d on it %d recieved not OK 2\n", e.rank, e.it);
							fixing_all(&e);
							continue;
						}
                	} 
           	 	}
            	
				flag = 1;
            	for (int i = 0; i < e.size; i++) 
				{
                	if (e.live_procs[i] && !eps_achieved[i]) 
					{
                    	flag = 0;
                    	break;
                	}
            	}
            	if (flag) 
				{
                	for (int i = 1; i < e.size; i++) 
					{
					   if (e.live_procs[i])
					   { 
						   rc = MPI_Isend(&flag, 1, MPI_INT, i,
								   EPS_TAG + itmax * err_shift, MPI_COMM_WORLD, &request);
						   if (rc != MPI_SUCCESS)
						   {
							   printf("Process %d on it %d recieved not OK 3\n", e.rank, e.it);
   							   fixing_all(&e);
   							   continue;
						   }
					   }
                	}
            	}
        	} 
			else 
			{
            	rc = MPI_Iprobe(0, EPS_TAG + itmax * err_shift, MPI_COMM_WORLD, &flag, &status);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not OK 4\n", e.rank, e.it);
					fixing_all(&e);
					continue;
				}	
				if (flag) 
				{
                	rc = MPI_Recv(&flag, 1, MPI_INT, 0, EPS_TAG + itmax * err_shift, MPI_COMM_WORLD, &status);
					if (rc != MPI_SUCCESS)
					{
						printf("Process %d on it %d recieved not OK 4\n", e.rank, e.it);
						fixing_all(&e);
						continue;
					}
                	should_calculate = 0;
           	 	}
        	}

    // get values
        	if (e.rank) 
			{ // process zero does not recieves data from other processes here
            //on first iteration left values could be counted by process itself
            	rc = MPI_Recv(matrix - N * N, N * N, MPI_DOUBLE, e.rank - 1, e.it + itmax * (err_shift - 1), 
						MPI_COMM_WORLD, &status);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not OK 5\n", e.rank, e.it);
					fixing_all(&e);
					continue;
				}
        	}
        	if (e.it && e.rank != e.size - 1) 
			{ // on first iteration next_values is already correct,
            //for last process next_values is ALWAYS correct
            	rc = MPI_Recv(matrix + N * N * (e.ind_last - e.ind_first), N * N, MPI_DOUBLE, 
						e.next_proc, e.it - 1 + itmax * (err_shift - 1), MPI_COMM_WORLD, &status);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not OK 6\n", e.rank, e.it);
					fixing_all(&e);
					continue;
				}
        	}
        // process values
        	eps = 0;
        	if (should_calculate) 
			{
        		for (ind1 = 0; ind1 < e.ind_last - e.ind_first; ind1++) 
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
			// backup counted messages
			char *filename = (char *)calloc(sizeof(char), 100);
			snprintf(filename, 100, "backup_%d_%d", e.rank, e.it);
			int fd = open(filename, O_CREAT, S_IRWXU);
		    if (!fd)
			{
				fprintf(stderr, "Process %d can not create and open file %s\n", e.rank, filename);
				MPI_Abort(MPI_COMM_WORLD, MPI_ERR_FILE);
			}	
			free(filename);
			write(fd, matrix, N * N * (e.ind_last - e.ind_first + 2) * sizeof(double));
			close(fd);
			// end of backup counted messages

			if (e.rank != e.size - 1)
			{
    	        rc = MPI_Isend(matrix + (e.ind_last - e.ind_first - 1) * N * N, N * N, MPI_DOUBLE, e.next_proc, e.it + itmax * (err_shift - 1), MPI_COMM_WORLD, &request);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not OK 7\n", e.rank, e.it);
					fixing_all(&e);
					continue;
				}
			}
        	if (e.rank && e.it != itmax - 1)
			{
            	rc = MPI_Isend(matrix, N * N, MPI_DOUBLE, e.prev_proc, e.it + itmax * (err_shift - 1), MPI_COMM_WORLD, &request);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not OK 8\n", e.rank, e.it);
					fixing_all(&e);
					continue;
				}				
			}
        	if (should_calculate)
			{
        		if (eps < maxeps)
				{
            // notify zero process about this
            		if (!e.rank)
					{
                		eps_achieved[0] = 1;
            		}
				   	else
					{
                		int achieved = 1;
                		rc = MPI_Send(&achieved, 1, MPI_INT, 0, EPS_TAG + itmax * err_shift, MPI_COMM_WORLD);
						if (rc != MPI_SUCCESS)
						{
							printf("Process %d on it %d recieved not OK 9\n", e.rank, e.it);
							fixing_all(&e);
							continue;
						}
            		}
         		}
        	}
			e.it++; // go to next iteration without problems
		}
//now we should sum all counted values and sent them to process 0
//process 0 should sum all recieved values and print them
    
	    double S = 0.;
   	 	for (ind1 = 0; ind1 < e.ind_last - e.ind_first; ind1++)
		{
        	for (ind2 = 1; ind2 < N - 1; ind2++)
			{
            	for (ind3 = 1; ind3 < N - 1; ind3++)
				{
                	S += *(matrix + ind1 * N * N + ind2 * N + ind3) * (ind3 + 1) * (ind2 + 1) * (ind1 + e.ind_first + 2) / (N * N * N);
            	}
        	}
    	}
    	if (e.rank)
		{
        	rc = MPI_Send(&S, 1, MPI_DOUBLE, 0, END_TAG + itmax * err_shift, MPI_COMM_WORLD);
			if (rc != MPI_SUCCESS)
			{
				printf("Process %d on it %d recieved not OK 10 NOT GOOD\n", e.rank, e.it);
				//fixing_all(&e); // out of cycle, all bad?
				//continue;
			}
    	}

    	if (!e.rank)
		{
        	double S_tmp;
        	for (ind1 = 1; ind1 < e.size; ind1++)
			{
				if (!e.live_procs[ind1])
				{
					continue;
				}
            	rc = MPI_Recv(&S_tmp, 1, MPI_DOUBLE, ind1, 
						END_TAG + itmax * err_shift, MPI_COMM_WORLD, &status);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not OK 11\n", e.rank, e.it);
					fixing_all(&e);
					continue;
				}
            	S += S_tmp;
        	}
    	}
    	if (!e.rank)
		{
			fprintf(stderr, "S = %lf\n", S);
    	}
    	free(matrix - N * N);
    }
    MPI_Finalize();
	return 0;
}
