#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

typedef struct particle_result_t {
	int is_right;
	int life_time;
} particle_result_t;

typedef struct result_t {
	double p;
	double average_life_time;
	double work_time;
} result_t;

typedef struct random_walk_ctx_t {
	int a;
	int b;
	int x;
	size_t N;
	double p;
	result_t res;
} random_walk_ctx_t;

unsigned int seed_for_threads[128];

int get_random_step(random_walk_ctx_t* ctx) {
	int thread = omp_get_thread_num();
	double random_value = (double) rand_r(&seed_for_threads[thread]) / RAND_MAX;
	if (random_value > ctx->p) {
		return -1;
	} else {
		return 1;
	}
}

int get_next(random_walk_ctx_t* ctx, int current) {
	return current + get_random_step(ctx);
}

particle_result_t process_particle(random_walk_ctx_t* ctx) {
	int position = ctx->x;
	int is_right = 0;
	int life_time = 0;
	while (1) {
		if (position <= ctx->a) {
			is_right = 0;
			break;
		} else if (position >= ctx->b) {
			is_right = 1;
			break;
		}
		position = get_next(ctx, position);
		++life_time;
	}
	particle_result_t res = {
		.is_right = is_right,
		.life_time = life_time,
	};
	return res;
}

void process_all_omp(random_walk_ctx_t* ctx) {
	int number_of_right = 0;
	int sum_life_time = 0;
	particle_result_t particle_result;

	double ts = omp_get_wtime();
	#pragma omp parallel for schedule(dynamic) reduction(+: number_of_right, sum_life_time)
	for (int i = 0; i < ctx->N; ++i) {
		particle_result = process_particle(ctx);
		number_of_right += particle_result.is_right;
		sum_life_time += particle_result.life_time;
	}
	double tf = omp_get_wtime();

	result_t res = {
		.p = (double) number_of_right / ctx->N,
		.average_life_time = (double) sum_life_time / ctx->N,
		.work_time = tf - ts,
	};
	ctx->res = res;
}

void process_all(random_walk_ctx_t* ctx) {
	int number_of_right = 0;
	int sum_life_time = 0;
	particle_result_t particle_result;

	double ts = omp_get_wtime();
	for (int i = 0; i < ctx->N; ++i) {
		particle_result = process_particle(ctx);
		number_of_right += particle_result.is_right;
		sum_life_time += particle_result.life_time;
	}
	double tf = omp_get_wtime();

	result_t res = {
		.p = (double) number_of_right / ctx->N,
		.average_life_time = (double) sum_life_time / ctx->N,
		.work_time = tf - ts,
	};
	ctx->res = res;
}

int main(int argc, char** argv) {
	if (argc < 7) {
		printf("There are not enough arguments\n");
		return 0;
	}

	int a = atoi(argv[1]);
	int b = atoi(argv[2]);
	int x = atoi(argv[3]);
	int N = atoi(argv[4]);
	double p = atof(argv[5]);
	int P = atoi(argv[6]);

	if (!(x > a && x < b)) {
		printf("x must be between a and b\n");
		return 0;
	}

	omp_set_num_threads(P);

	srand(time(NULL));

	for (int i = 0; i < P; i++) {
		seed_for_threads[i] = rand();
	}

	random_walk_ctx_t ctx = {
		.a = a,
		.b = b,
		.x = x,
		.N = N,
		.p = p,
	};

	if (P < 2) {
		process_all(&ctx);
	} else {
		process_all_omp(&ctx);
	}

	FILE* file = fopen("stats.txt", "w");
	if (file == NULL) {
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(file, "%.2f %.1f %.5fs %d %d %d %d %.2f %d\n", ctx.res.p, ctx.res.average_life_time, 
		ctx.res.work_time, a, b, x, N, p, P);

	fclose(file);

	return 0;
}
