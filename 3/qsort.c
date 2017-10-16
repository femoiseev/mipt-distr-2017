#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

int compare (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Incorrect number of arguments\n");
    }

    int n = atoi(argv[1]);
    int* ar = (int*) malloc(n * sizeof(int));
    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        ar[i] = rand();
    }

    double start = omp_get_wtime();
    qsort(ar, n, sizeof(int), compare);
    double end = omp_get_wtime();
    printf("%.5fs\n", end - start);

    free(ar);

    return EXIT_SUCCESS;
}