CC = mpicc
CFLAGS = -std=c99
LDFLAGS = -fopenmp

all:
	$(CC) $(CFLAGS) $(LDFLAGS) -o rw rw.c

.PHONY: run
run: all
	mpirun -np 4 ./rw 1 2 2 100 10 0.25 0.25 0.25 0.25

batch: all
	sbatch -n 4 ./wrapper ./rw 1 2 2 100 10 0.25 0.25 0.25 0.25
