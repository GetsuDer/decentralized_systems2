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
	int *alive_procs; // array of currently live processes
	int err_shift; // value for tag shifting
	int reread_it; // value to reread iteration
};

// rearrange_buff
// [0]: dead_proc ind
// [1]: it ind
// [2]: is_slave
void
fixing_all(struct env *env, int err, int *rearrange_buff)
{
	printf("In fixing_all proc %d/%d on it %d\n", env->rank, env->size, env->it);
	int dead_proc = rearrange_buff[0];
	int is_slave = rearrange_buff[2];
	int buff[3];
	MPI_Status status;
	if (!env->alive_procs[dead_proc])
	{
// повторное сообщение
		printf("Proc %d recieved second message about death of %d", env->rank, dead_proc);
		return;
	}
	env->alive_procs[dead_proc] = 0;
	// пересчитать значения
	int live_number = 1;
	int prev_live = -1;
	int next_live = -1;
	for (int i = env->rank - 1; i >= 0; i--)
	{
		if (env->alive_procs[i])
		{
			prev_live = i;
			break;
		}
	}
	if (prev_live == -1)
	{
		for (int i = env->size - 1; i > env->rank; i--)
		{
			if (env->alive_procs[i])
			{
				prev_live = i;
				break;
			}
		}
	}	
	for (int i = env->rank + 1; i < env->size; i++)
	{
		if (env->alive_procs[i])
		{
			next_live = i;
			break;
		}
	}
	if (next_live == -1)
	{
		for (int i = 0; i < env->rank; i++)
		{
			if (env->alive_procs[i])
			{
				next_live = i;
				break;
			}
		}
	}
	int new_rank = -1;
	for (int i = 0; i < env->size; i++)
	{
		if (env->alive_procs[i])
		{
			if (i <= env->rank)
			{
				new_rank++;
			}
			live_number++;
		}
	}
	int i_am_prev = (env->prev_proc == prev_live);
	if (!is_slave)
	{
		if (i_am_prev) {
			MPI_Recv(buff, 3, MPI_INT, MPI_ANY_SOURCE, RESTORE_TAG + itmax * env->err_shift, MPI_COMM_WORLD, &status);
			env->reread_it = buff[1];
			is_slave = 1;
		}
	}
  	int part_size = (N - 2) / live_number;
	env->ind_last = new_rank * part_size; 
	env->ind_first = env->ind_last + part_size;
	if (next_live < env->rank)
	{ // last proc
		env->ind_last = N - 2;
	}
	env->prev_proc = prev_live;
	env->next_proc = next_live;

	// restore iteration and continue
	if (is_slave)
	{
		env->reread_it = rearrange_buff[1];
		env->err_shift++;
		return;
	}
	else
	{
		printf("Proc %d starting restoring.\n", env->rank);

		int max_good_it = env->it - 1; // for this proc previous it is backuped for sure
		for (int i = 0; i < env->size; i++)
		{
			if (env->alive_procs[i] && (i != env->rank))
			{
				for (; max_good_it > 0; max_good_it--)
				{
					char *filename = (char *)calloc(100, sizeof(char));
					snprintf(filename, 100, "backup_%d_%d", i, max_good_it);
					if (!access(filename, F_OK)) // file exists
					{
						free(filename);
						snprintf(filename, 100, "backup_%d_%d", i, max_good_it + 1);
						if (access(filename, F_OK)) // file not exists
						{
							max_good_it--;
						}
						break;
					}
				}
			}
		}
		printf("Founded correct it: %d\n", max_good_it);
		buff[2] = 1;
		buff[0] = dead_proc;
		buff[1] = max_good_it;
		for (int i = 0; i < env->size; i++)
		{
			if (env->alive_procs[i] && (i != env->rank))
			{
				printf("%d sending restore message to %d\n", env->rank, i);
				MPI_Send(buff, 3, MPI_INT, i, 
						RESTORE_TAG + itmax * env->err_shift, MPI_COMM_WORLD);
			}
		}
		env->err_shift++;
		env->reread_it = max_good_it;
		printf("Master, continuing normal work\n");
		return;
	}
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
	e.err_shift = 1;
	e.reread_it = -1;
	int rearrange_buff[3];
	MPI_Request request;

    MPI_Comm_size(MPI_COMM_WORLD, &e.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &e.rank);

	if (!(bad_rank % e.size))
	{
		bad_rank++;
	}

	if (!e.rank) 
		printf("bad_it %d bad_rank %d\n", bad_it, bad_rank % e.size);

    if (e.size > N - 2) 
	{
        e.size = N - 2;
    }

    if (e.rank < e.size) 
	{
		e.alive_procs = (int *)calloc(e.size, sizeof(int));
		for (int i = 0; i < e.size; i++)
		{
			e.alive_procs[i] = 1;
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
		if (e.ind_last != N - 2) 
		{ //for last process next_values is always zero submatrix
        	for (ind2 = 0; ind2 <= N - 1; ind2++) 
			{
            	for (ind3 = 0; ind3 <= N - 1; ind3++) 
				{
                	if (ind2 == 0 || ind3 == 0 || ind3 == N - 1 || ind2 == N - 1 ) {
                    	*(matrix + N * N * (e.ind_last - e.ind_first) + ind2 * N + ind3) = 0;
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
			int restoring = 0;
			if ((e.it == bad_it) && e.rank == (bad_rank % e.size))
			{
				// raise error
				printf("On iteration %d process %d/%d decided to die.\n", e.it, e.rank, e.size);
				exit(1);
			}

			if (e.reread_it != -1)
			{
				// return to one of ther previous iterations
				printf("Process %d rereading matrix for it %d\n", e.rank, e.reread_it);
				char *filename = (char *)calloc(sizeof(char), 100);
				snprintf(filename, 100, "backup_%d_%d", e.rank, e.reread_it);
				int fd = open(filename, O_RDONLY);
				if (!fd)
				{
					printf("Proc %d can not open backup file for it %d\n", e.rank, e.reread_it);
					MPI_Abort(MPI_COMM_WORLD, MPI_ERR_FILE);
				}
				free(filename);
				read(fd, matrix, N * N * (e.ind_last - e.ind_first + 2) * sizeof(double));
				printf("Process %d restored and continuing work on it %d\n", e.rank, e.reread_it);
				e.it = e.reread_it;
				e.reread_it = -1;
				restoring = 1;
			}
	
        	int flag;
        	if (!e.rank) 
			{
// check, if have messages about eps achievement from other processes
				int all_good = 1;
            	for (int process_id = 1; all_good && (process_id < e.size); process_id++) 
				{
					if (!e.alive_procs[process_id])
					{
						// this process is already dead, forget about it
						continue;
					}
                	flag = 0;
                	rc = MPI_Iprobe(process_id, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
					if (rc != MPI_SUCCESS)
					{
						printf("Process %d on it %d recieved not OK from source %d\n", e.rank, e.it, status.MPI_SOURCE);
						rearrange_buff[2] = 0;
						rearrange_buff[0] = process_id;
						fixing_all(&e, rc, rearrange_buff);
						continue; // go to next iteration
					}
					
					if (flag)
					{
						if (status.MPI_TAG == (EPS_TAG + itmax * e.err_shift))
						{

							eps_achieved[process_id] = 1;
		                    rc = MPI_Recv(&flag, 1, MPI_INT, process_id, 
									EPS_TAG + itmax * e.err_shift, MPI_COMM_WORLD, &status); 
								// take this
                	    	//message from queue
							if (rc != MPI_SUCCESS)
							{
								printf("Process %d on it %d recieved not OK from %d\n", e.rank, e.it, process_id);
								rearrange_buff[2] = 0;
								rearrange_buff[0] = process_id;
								fixing_all(&e, rc, rearrange_buff);
								all_good = 0;
								break;
							}
						}
						else if (status.MPI_TAG == (RESTORE_TAG + itmax * e.err_shift))
						{
		                    rc = MPI_Recv(&rearrange_buff, 3, MPI_INT, process_id, 
									RESTORE_TAG + itmax * e.err_shift, MPI_COMM_WORLD, &status); 
							printf("Process 0 recieve restore tag from process %d\n", process_id);
							rearrange_buff[2] = 1;
							rearrange_buff[0] = process_id;
							fixing_all(&e, rc, rearrange_buff);
							all_good = 0; //to go to next iteration
						}
                	} 
           	 	}
                if (!all_good)
				{
					continue;
				}

				flag = 1;
            	for (int i = 0; i < e.size; i++) 
				{
                	if (e.alive_procs[i] && !eps_achieved[i]) 
					{
                    	flag = 0;
                    	break;
                	}
            	}
            	if (flag) 
				{
                	for (int i = 1; i < e.size; i++) 
					{
					   if (e.alive_procs[i])
					   { 
						   rc = MPI_Isend(&flag, 1, MPI_INT, i,
								   EPS_TAG + itmax * e.err_shift, MPI_COMM_WORLD, &request);
						   if (rc != MPI_SUCCESS)
						   {
							   printf("Process %d on it %d recieved not OK from proc %d\n", e.rank, e.it, i);
   							   rearrange_buff[0] = i;
							   rearrange_buff[2] = 0;
							   fixing_all(&e, rc, rearrange_buff);
   							   continue;
						   }
					   }
                	}
            	}
        	} 
			else 
			{
            	rc = MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not from process 0\n", e.rank, e.it);
					MPI_Abort(MPI_COMM_WORLD, MPI_ERR_FILE);
				}	
				if (flag) 
				{
					if (status.MPI_TAG == (EPS_TAG + itmax * e.err_shift))
					{
						rc = MPI_Recv(&flag, 1, MPI_INT, 0, 
								EPS_TAG + itmax * e.err_shift, MPI_COMM_WORLD, &status);
						if (rc != MPI_SUCCESS)
						{
							printf("Process %d on it %d recieved not OK from process 0\n", e.rank, e.it);
							MPI_Abort(MPI_COMM_WORLD, MPI_ERR_FILE);
						}
                		should_calculate = 0;
					}
					else if (status.MPI_TAG == (RESTORE_TAG + itmax * e.err_shift))
					{
						rc = MPI_Recv(rearrange_buff, 3, MPI_INT, status.MPI_SOURCE, 
							   	RESTORE_TAG + itmax * e.err_shift, MPI_COMM_WORLD, &status);
						printf("Process %d revieved restore tag from proc %d\n", e.rank, status.MPI_SOURCE);
						rearrange_buff[2] = 1;
						rearrange_buff[0] = status.MPI_SOURCE;
						fixing_all(&e, rc, rearrange_buff);
					}	
           	 	}
        	}

    // get values
        	if (e.rank) 
			{ // process zero does not recieves data from other processes here
            //on first iteration left values could be counted by process itself
				int go_to_next = 0;
				printf("proc %d waiting for left data from proc %d on it %d with tag %d\n", 
						e.rank, e.prev_proc, e.it, e.it + itmax * (e.err_shift - 1));
				while (1) {
					// если есть сообщение о восстановлении, берем его и идем на следующий цикл 
            		rc = MPI_Iprobe(MPI_ANY_SOURCE, RESTORE_TAG + itmax * e.err_shift, MPI_COMM_WORLD, &flag, &status);
					if (flag) {
						MPI_Recv(&rearrange_buff, 3, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
						fixing_all(&e, rc, rearrange_buff);
						go_to_next = 1;
						break;
					}
					// иначе должно быть сообщение от другого процесса, пока его нет, повторяем ожидание
					rc = MPI_Iprobe(e.prev_proc, e.it + itmax * (e.err_shift - 1), MPI_COMM_WORLD, &flag, &status);	
					if (rc != MPI_SUCCESS)
					{
						printf("Process %d on it %d recieved not OK from proc %d\n", e.rank, e.it, e.prev_proc);
						rearrange_buff[2] = 0;
						rearrange_buff[0] = e.prev_proc; 
						fixing_all(&e, rc, rearrange_buff);
						go_to_next = 1;
						break;
					}
					if (flag) {
						rc = MPI_Recv(matrix - N * N, N * N, MPI_DOUBLE, 
							e.prev_proc, e.it + itmax * (e.err_shift - 1), MPI_COMM_WORLD, &status);
						if (rc != MPI_SUCCESS)
						{
							printf("Process %d on it %d recieved not OK from proc %d\n", e.rank, e.it, e.prev_proc);
							rearrange_buff[2] = 0;
							rearrange_buff[0] = e.prev_proc; 
							fixing_all(&e, rc, rearrange_buff);
							go_to_next = 1;
							break;
						}
						break;
					}
				}
        		if (go_to_next)
				{
					continue;
				}
			}
			printf("proc %d recieved left data from proc %d on it %d with tag %d\n", 
					e.rank, e.prev_proc, e.it, e.it + itmax * (e.err_shift - 1));
        	if ((e.it && !restoring) && e.rank != e.size - 1) 
			{ // on first iteration next_values is already correct,
            //for last process next_values is ALWAYS correct
				int go_to_next = 0;
				printf("proc %d waiting for right data from proc %d on it %d with tag %d\n", 
						e.rank, e.next_proc, e.it - 1, e.it - 1 + itmax * (e.err_shift - 1));
				while (1) {
					// если есть сообщение о восстановлении, берем его и идем на следующий цикл 
            		rc = MPI_Iprobe(MPI_ANY_SOURCE, RESTORE_TAG + itmax * e.err_shift, MPI_COMM_WORLD, &flag, &status);
					if (flag) {
						MPI_Recv(&rearrange_buff, 3, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
						fixing_all(&e, rc, rearrange_buff);
						go_to_next = 1;
						break;
					}
					// иначе должно быть сообщение от другого процесса, пока его нет, повторяем ожидание
					rc = MPI_Iprobe(e.next_proc, e.it - 1 + itmax * (e.err_shift - 1), MPI_COMM_WORLD, &flag, &status);	
					if (rc != MPI_SUCCESS)
					{
						printf("Process %d on it %d recieved not OK from proc %d 496\n", e.rank, e.it, e.next_proc);
						rearrange_buff[2] = 0;
						rearrange_buff[0] = e.next_proc; 
						fixing_all(&e, rc, rearrange_buff);
						go_to_next = 1;
						break;
					}
					if (flag) {
            			rc = MPI_Recv(matrix + N * N * (e.ind_last - e.ind_first), N * N, MPI_DOUBLE, 
							e.next_proc, e.it - 1 + itmax * (e.err_shift - 1), MPI_COMM_WORLD, &status);
						if (rc != MPI_SUCCESS)
						{
							printf("Process %d on it %d recieved not OK from proc %d 508\n", e.rank, e.it, e.next_proc);
							rearrange_buff[2] = 0;
							rearrange_buff[0] = e.next_proc; 
							fixing_all(&e, rc, rearrange_buff);
							go_to_next = 1;
							break;
						}
						break;
					}
				}
				if (go_to_next)
				{
					continue;
				}
        	}
			printf("proc %d recieved right data from proc %d on it %d with tag %d\n", 
					e.rank, e.prev_proc, e.it - 1, e.it - 1 + itmax * (e.err_shift - 1));
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
			printf("Process %d ended it %d\n", e.rank, e.it);
			// end of backup counted messages

			if (e.rank != e.size - 1)
			{
				printf("Proc %d sending message to %d with tag %d\n", e.rank, e.next_proc, e.it + itmax * (e.err_shift - 1));
    	        rc = MPI_Isend(matrix + (e.ind_last - e.ind_first - 1) * N * N, N * N, MPI_DOUBLE, e.next_proc, e.it + itmax * (e.err_shift - 1), MPI_COMM_WORLD, &request);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not OK from proc %d\n", e.rank, e.it, e.next_proc);
					rearrange_buff[0] = e.next_proc;
					rearrange_buff[2] = 0;
					fixing_all(&e, rc, rearrange_buff);
					continue;
				}
			}
        	if (e.rank && e.it != itmax - 1)
			{
				printf("Proc %d sending message to %d with tag %d\n", e.rank, e.prev_proc, e.it + itmax * (e.err_shift - 1));
            	rc = MPI_Isend(matrix, N * N, MPI_DOUBLE, e.prev_proc, e.it + itmax * (e.err_shift - 1), MPI_COMM_WORLD, &request);
				if (rc != MPI_SUCCESS)
				{
					printf("Process %d on it %d recieved not OK from proc %d\n", e.rank, e.it, e.prev_proc);
					rearrange_buff[0] = e.prev_proc;
					rearrange_buff[2] = 0;
					fixing_all(&e, rc, rearrange_buff);
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
                		rc = MPI_Send(&achieved, 1, MPI_INT, 0, EPS_TAG + itmax * e.err_shift, MPI_COMM_WORLD);
						if (rc != MPI_SUCCESS)
						{
							printf("Process %d on it %d recieved not OK from proc 0 (this is strange 3)\n", e.rank, e.it);
							MPI_Abort(MPI_COMM_WORLD, MPI_ERR_FILE);
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
        	rc = MPI_Send(&S, 1, MPI_DOUBLE, 0, END_TAG + itmax * e.err_shift, MPI_COMM_WORLD);
    	}

    	if (!e.rank)
		{
        	double S_tmp;
        	for (ind1 = 1; ind1 < e.size; ind1++)
			{
				if (!e.alive_procs[ind1])
				{
					continue;
				}
            	rc = MPI_Recv(&S_tmp, 1, MPI_DOUBLE, ind1, 
						END_TAG + itmax * e.err_shift, MPI_COMM_WORLD, &status);
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
