CC1=gcc -w -O3
CC2=nvcc -w -O3 

EXECUTABLES = serial_meanshift cuda_meanshift_gm cuda_meanshift_sm

all: $(EXECUTABLES)

serial_meanshift: serial_meanshfit.c
	$(CC1) $< -o $@	

cuda_meanshift_gm: cuda_meanshift_gm.cu
	$(CC2) $< -o $@

cuda_meanshift_sm: cuda_meanshift_sm.cu
	$(CC2) $< -o $@

