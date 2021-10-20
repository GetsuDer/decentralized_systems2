SRCDIR=src/
MPICC=mpicc
NPROC=10
.PHONY: all clean mpi

all: mpi
	mpiexec -n $(NPROC) ./mpi

mpi: $(SRCDIR)mpi_code.c
	$(MPICC) $(SRCDIR)mpi_code.c -o mpi

clean:
	rm -f mpi
