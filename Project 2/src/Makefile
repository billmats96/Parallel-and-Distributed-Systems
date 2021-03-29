CC1=gcc -O3
CC2=mpicc -O3
CC3=mpicc -O3 -fopenmp

EXECUTABLES = mpi_blocking mpi_nonblocking mpi_threaded_blocking mpi_threaded_nonblocking serialKNN

all: $(EXECUTABLES)

serialKNN: serialKNN.c
	$(CC1) $< -o $@	

mpi_blocking: mpi_blocking.c
	$(CC2) $< -o $@

mpi_nonblocking: mpi_nonblocking.c
	$(CC2) $< -o $@

mpi_threaded_blocking: mpi_threaded_blocking.c
	$(CC3) $< -o $@

mpi_threaded_nonblocking: mpi_threaded_nonblocking.c
	$(CC3) $< -o $@

