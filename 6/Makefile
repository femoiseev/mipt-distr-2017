build:
	mpicc rw.c -o run -fopenmp

build_single:
	mpicc rw_single.c -o single -fopenmp

work:
	mpirun -np 4 ./run 1 2 2 100 10 0.25 0.25 0.25 0.25

work_single:
	mpirun -np 4 ./single 1 2 2 100 10 0.25 0.25 0.25 0.25