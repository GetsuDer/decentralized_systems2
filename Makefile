SRCDIR=src/
#SRCDIR=.
MPICC=mpicc
NPROC=10
.PHONY: all clean mpi mpi_upgraded

all: mpi mpi_upgraded

run_mpi: mpi
	mpiexec -n $(NPROC) ./mpi

run_mpi_upgraded: mpi_upgraded
	rm -f backup_*
	mpiexec -n $(NPROC) --disable-auto-cleanup ./mpi_upgraded

mpi_upgraded:
	$(MPICC) $(SRCDIR)mpi_upgraded.c -o mpi_upgraded

mpi: $(SRCDIR)mpi_code.c
	$(MPICC) $(SRCDIR)mpi_code.c -o mpi

clean:
	rm -f mpi mpi_upgraded backup_*
