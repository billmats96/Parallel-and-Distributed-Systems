CC1=gcc -w -O3 -lm
CC2=nvcc -w -O3 

EXECUTABLES = serialnlm serialnlm_adaptive nlm adaptive 

all: $(EXECUTABLES)

serialnlm: serialnlm.c
	$(CC1) $< -o $@	

serialnlm_adaptive: serialnlm_adaptive.c
	$(CC1) $< -o $@

nlm: nlm.cu
	$(CC2) $< -o $@

adaptive: adaptive.cu
	$(CC2) $< -o $@

