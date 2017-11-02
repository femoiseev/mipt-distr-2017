#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#define MASTER 0
#define UP 0
#define RIGHT 1
#define DOWN 2
#define LEFT 3
#define CONN_TIME 100

typedef struct ctx_t {
    int l;
    int a;
    int b;
    int n;
    int N;
    double pl;
    double pr;
    double pu;
    double pd;
} ctx_t;

typedef struct particle_t {
    int x;
    int y;
    int n;
} particle_t;

void swap_particle_array(particle_t** x, particle_t** y) {
    particle_t* tmp = *x;
    *x = *y;
    *y = tmp;
}

void swap_int(int* x, int* y) {
    int tmp = *x;
    *x = *y;
    *y = tmp;
}


int get_dir(double left, double up, double right, double down) {
    if (left >= up && left >= right && left >= down) {
        return LEFT;
    } else if (up >= left && up >= right && up >= down) {
        return UP;
    } else if (right >= left && right >= up && right >= down) {
        return RIGHT;
    } else {
        return DOWN;
    }
}

int step(ctx_t* ctx) {
    double left = rand() * ctx->pl;
    double right = rand() * ctx->pr;
    double up = rand() * ctx->pu;
    double down = rand() * ctx->pd;

    return get_dir(left, up, right, down);
}

void insert(particle_t x, particle_t** ar, int* size, int* max_size) {
    if (*size >= *max_size) {
        *max_size *= 2;
        *ar = (particle_t* )realloc(*ar, (*max_size) * sizeof(particle_t));
    }
    (*ar)[*size] = x;
    (*size)++;
}

void delete(int index, particle_t** ar, int* size) {
    (*ar)[index] = (*ar)[(*size)-1];
    (*size)--;
}

void process(ctx_t* ctx, int rank, int comm_size) {
    int a = rank % ctx->a;
    int b = rank / ctx->a;
    int left_rank, right_rank, up_rank, down_rank;
    if (rank % ctx->a != 0) {
        left_rank = rank - 1;
    } else {
        left_rank = rank + ctx->a - 1;
    }

    if (rank % ctx->a != ctx->a - 1) {
        right_rank = rank + 1;
    } else {
        right_rank = rank - (ctx->a - 1);
    }

    if (rank / ctx->a != 0) {
        up_rank = rank - ctx->a;
    } else {
        up_rank = comm_size - ctx->a + rank;
    }

    if (rank / ctx->a != ctx->b-1) {
        down_rank = rank + ctx->a;
    } else {
        down_rank = rank % ctx->a;
    }

    int size = ctx->N;
    int max_size = 2 * size;
    particle_t* particles = (particle_t*) malloc(max_size * sizeof(particle_t));

    int left_size = 0;
    int right_size = 0;
    int up_size = 0;
    int down_size = 0;
    int left_max_size = size;
    int right_max_size = size;
    int up_max_size = size;
    int down_max_size = size;
    particle_t* to_left = (particle_t*) malloc(left_max_size * sizeof(particle_t));
    particle_t* to_right = (particle_t*) malloc(right_max_size * sizeof(particle_t));
    particle_t* to_up = (particle_t*) malloc(up_max_size * sizeof(particle_t));
    particle_t* to_down = (particle_t*) malloc(down_max_size * sizeof(particle_t));

    int fin_size = 0;
    int fin_max_size = size;
    particle_t* finished = (particle_t*) malloc(fin_max_size * sizeof(particle_t));
    
    omp_lock_t lock;
    int is_running = 1;
    omp_init_lock(&lock);
    omp_set_lock(&lock);

    #pragma omp parallel sections default(shared)
    {
        #pragma omp section
        {
            while(is_running) {
                omp_set_lock(&lock);
                //printf("1: %d\n", rank);
                for (int j = 0; j < CONN_TIME; j++) {
                    int i = 0;
                    while(i < size) {
                    
                    //for (int i = 0; i < size; i++) {
                        if (particles[i].n == 0) {
                            //printf("%d: prev fin size: %d, prev size: %d, i: %d\n", rank, fin_size, size, i);
                            insert(particles[i], &finished, &fin_size, &fin_max_size);
                            delete(i, &particles, &size);
                            //printf("%d: fin size: %d, size: %d\n", rank, fin_size, size);
                            continue;
                        }
                        //printf("rank: %d, i: %d, n: %d\n", rank, i, particles[i].n);
                        particles[i].n--;
                        int dir = step(ctx);
                        switch(dir) {
                            case LEFT:
                                particles[i].x--;
                                if (particles[i].x < 0) {
                                    particles[i].x = ctx->l + particles[i].x;
                                    insert(particles[i], &to_left, &left_size, &left_max_size);
                                    delete(i, &particles, &size);
                                    i--;
                                }
                                break;
                            case RIGHT:
                                particles[i].x++;
                                if (particles[i].x >= ctx->l) {
                                    particles[i].x = ctx->l - particles[i].x;
                                    insert(particles[i], &to_right, &right_size, &right_max_size);
                                    delete(i, &particles, &size);
                                    i--;
                                }
                                break;
                            case UP:
                                particles[i].y--;
                                if (particles[i].y < 0) {
                                    particles[i].y = ctx->l + particles[i].y;
                                    insert(particles[i], &to_up, &up_size, &up_max_size);
                                    delete(i, &particles, &size);
                                    i--;
                                }
                                break;
                            case DOWN:
                                particles[i].y++;
                                if (particles[i].y >= ctx->l) {
                                    particles[i].y = ctx->l - particles[i].y;
                                    insert(particles[i], &to_down, &down_size, &down_max_size);
                                    delete(i, &particles, &size);
                                    i--;
                                }
                                break;
                            default:
                                break;
                        }
                        i++;
                        //printf("rank: %d, i: %d, n: %d\n", rank, i, particles[i].n);
                    }
                }
                omp_unset_lock(&lock);
                //#pragma omp barrier
            }
        }

        #pragma omp section
        {
            int* seeds;
            int seed;
            if (rank == MASTER) {
                srand(time(NULL));
                seeds = (int*) malloc (comm_size * sizeof(int));
                for (int i = 0; i < comm_size; i++) {
                    seeds[i] = rand();
                }
            }
            MPI_Scatter(seeds, 1, MPI_INT, &seed, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
            srand(seed);

            for (int i = 0; i < size; i++) {
                int x = rand() % ctx->l;
                int y = rand() % ctx->l;
                particle_t tmp_particle = {
                    .x = x,
                    .y = y,
                    .n = ctx->n,
                };
                particles[i] = tmp_particle;
            }

            int tmp_left_size = left_size;
            int tmp_right_size = right_size;
            int tmp_up_size = up_size;
            int tmp_down_size = down_size;
            int tmp_left_max_size = left_size;
            int tmp_right_max_size = right_size;
            int tmp_up_max_size = up_size;
            int tmp_down_max_size = down_size;

            particle_t* tmp_left = malloc(tmp_left_max_size * sizeof(int));
            particle_t* tmp_right = malloc(tmp_right_max_size * sizeof(int));
            particle_t* tmp_up = malloc(tmp_up_max_size * sizeof(int));
            particle_t* tmp_down = malloc(tmp_down_max_size * sizeof(int));
            int tmp_finished_size = fin_size;

            omp_unset_lock(&lock);

            while (is_running) {
                omp_set_lock(&lock);
                swap_int(&tmp_left_size, &left_size);
                swap_int(&tmp_right_size, &right_size);
                swap_int(&tmp_up_size, &up_size);
                swap_int(&tmp_down_size, &down_size);
                swap_int(&tmp_left_max_size, &left_max_size);
                swap_int(&tmp_right_max_size, &right_max_size);
                swap_int(&tmp_up_max_size, &up_max_size);
                swap_int(&tmp_down_max_size, &down_max_size);
                swap_particle_array(&tmp_left, &to_left);
                swap_particle_array(&tmp_right, &to_right);
                swap_particle_array(&tmp_up, &to_up);
                swap_particle_array(&tmp_down, &to_down);
                tmp_finished_size = fin_size;

                //omp_unset_lock(&lock);

                MPI_Request requests[8];

                int from_left_size, from_right_size, from_up_size, from_down_size;
                MPI_Issend(&tmp_left_size, 1, MPI_INT, left_rank, LEFT, MPI_COMM_WORLD, requests);
                MPI_Issend(&tmp_right_size, 1, MPI_INT, right_rank, RIGHT, MPI_COMM_WORLD, requests + 1);
                MPI_Issend(&tmp_up_size, 1, MPI_INT, up_rank, UP, MPI_COMM_WORLD, requests + 2);
                MPI_Issend(&tmp_down_size, 1, MPI_INT, down_rank, DOWN, MPI_COMM_WORLD, requests + 3);

                MPI_Irecv(&from_left_size, 1, MPI_INT, left_rank, RIGHT, MPI_COMM_WORLD, requests + 4);
                MPI_Irecv(&from_right_size, 1, MPI_INT, right_rank, LEFT, MPI_COMM_WORLD, requests + 5);
                MPI_Irecv(&from_up_size, 1, MPI_INT, up_rank, DOWN, MPI_COMM_WORLD, requests + 6);
                MPI_Irecv(&from_down_size, 1, MPI_INT, down_rank, UP, MPI_COMM_WORLD, requests + 7);

                MPI_Waitall(8, requests, MPI_STATUS_IGNORE);

                particle_t* from_left = (particle_t*) malloc(from_left_size * sizeof(particle_t));
                particle_t* from_right = (particle_t*) malloc(from_right_size * sizeof(particle_t));
                particle_t* from_up = (particle_t*) malloc(from_up_size * sizeof(particle_t));
                particle_t* from_down = (particle_t*) malloc(from_down_size * sizeof(particle_t));

                MPI_Issend(tmp_left, tmp_left_size * sizeof(particle_t), MPI_BYTE, left_rank, LEFT, MPI_COMM_WORLD, requests);
                MPI_Issend(tmp_right, tmp_right_size * sizeof(particle_t), MPI_BYTE, right_rank, RIGHT, MPI_COMM_WORLD, requests + 1);
                MPI_Issend(tmp_up, tmp_up_size * sizeof(particle_t), MPI_BYTE, up_rank, UP, MPI_COMM_WORLD, requests + 2);
                MPI_Issend(tmp_down, tmp_down_size * sizeof(particle_t), MPI_BYTE, down_rank, DOWN, MPI_COMM_WORLD, requests + 3);

                MPI_Irecv(from_left, from_left_size * sizeof(particle_t), MPI_BYTE, left_rank, RIGHT, MPI_COMM_WORLD, requests + 4);
                MPI_Irecv(from_right, from_right_size * sizeof(particle_t), MPI_BYTE, right_rank, LEFT, MPI_COMM_WORLD, requests + 5);
                MPI_Irecv(from_up, from_up_size * sizeof(particle_t), MPI_BYTE, up_rank, DOWN, MPI_COMM_WORLD, requests + 6);
                MPI_Irecv(from_down, from_down_size * sizeof(particle_t), MPI_BYTE, down_rank, UP, MPI_COMM_WORLD, requests + 7);

                MPI_Waitall(8, requests, MPI_STATUS_IGNORE);

                int all_finished[comm_size];
                MPI_Gather(&tmp_finished_size, 1, MPI_INT, all_finished, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

                //omp_set_lock(&lock);

                for (int i = 0; i < from_left_size; i++) {
                    insert(from_left[i], &particles, &size, &max_size);
                }

                for (int i = 0; i < from_right_size; i++) {
                    insert(from_right[i], &particles, &size, &max_size);
                }

                for (int i = 0; i < from_up_size; i++) {
                    insert(from_up[i], &particles, &size, &max_size);
                }

                for (int i = 0; i < from_down_size; i++) {
                    insert(from_down[i], &particles, &size, &max_size);
                }

                int is_actives[comm_size];

                if (rank == MASTER) {

                    int sum_finished = 0;
                    for (int i = 0; i < comm_size; i++) {
                        //printf("%d ", all_finished[i]);
                        sum_finished += all_finished[i];
                    }

                    if (sum_finished == comm_size * ctx->N) {
                        for (int i = 0; i < comm_size; i++) {
                            is_actives[i] = 0;
                        }
                    } else {
                        for (int i = 0; i < comm_size; i++) {
                            is_actives[i] = 1;
                        }
                    }

                    for (int i = 0; i < comm_size; i++) {
                        printf("%d ", is_actives[i]);
                    }
                    printf("\n");
                }

                //printf("%d: %d\n", rank, is_running);
                MPI_Scatter(is_actives, 1, MPI_INT, &is_running, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

                omp_unset_lock(&lock);

                //#pragma omp_barrier

                free(from_left);
                free(from_right);
                free(from_up);
                free(from_down);
            }

            printf("%d: %d\n", rank, tmp_finished_size);

            free(tmp_left);
            free(tmp_right);
            free(tmp_up);
            free(tmp_down);
        }
    }

    free(particles);
    free(to_left);
    free(to_right);
    free(to_up);
    free(to_down);
    omp_destroy_lock(&lock);
}

int main (int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 10) {
        if (rank == MASTER) {
            printf("Incorrect number of arguments\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    ctx_t ctx = {
        .l = atoi(argv[1]),
        .a = atoi(argv[2]),
        .b = atoi(argv[3]),
        .n = atoi(argv[4]),
        .N = atoi(argv[5]),
        .pl = atof(argv[6]),
        .pr = atof(argv[7]),
        .pu = atof(argv[8]),
        .pd = atof(argv[9]),
    };

    if (ctx.a * ctx.b != size) {
        if (rank == MASTER) {
            printf("Incorrect number of processes\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    omp_set_num_threads(2);

    process(&ctx, rank, size);

    MPI_Finalize();
    return 0;
}