CC = mpicc
CFLAGS = -std=c99
LDFLAGS = -fopenmp

all: rw

rw: rw.o

rw_single: rw_single.o

.PHONY: run
run: rw
        mpirun -np 4 ./rw 1 2 2 100 10 0.25 0.25 0.25 0.25

batch: rw
        sbatch -n 4 ./wrapper ./rw 1 2 2 100 10 0.25 0.25 0.25 0.25