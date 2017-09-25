#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Incorrect number of arguments");
        return EXIT_FAILURE;
    }

    int n =  atoi(argv[1]);
    srand(time(NULL));

    FILE* input = fopen("input.txt", "w");

    for (int i = 0; i < n; i++) {
        int val = rand();
        fprintf(input, "%d ", val);
    }
    fprintf(input, "\n");

    fclose(input);

    return EXIT_SUCCESS;
}