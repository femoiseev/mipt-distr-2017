#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <omp.h>
#include <string.h>

typedef struct sort_ctx_t {
    int n;
    int m;
    int* array;
    int* tmp;
    int threads_num;
    pthread_t* threads;
} sort_ctx_t;

typedef struct merge_ctx_t {
    int left_size;
    int right_size;
    int* left;
    int* right;
    int* target;
    int* tmp;
    int m;
    int threads_num;
    pthread_t* threads;
} merge_ctx_t;

int compare (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}

void merge(merge_ctx_t* ctx) {
    int left_index = 0;
    int right_index = 0;
    int index = 0;

    while (left_index < ctx->left_size && right_index < ctx->right_size) {
        if (ctx->left[left_index] < ctx->right[right_index]) {
            ctx->tmp[index] = ctx->left[left_index];
            ++left_index;
        } else {
            ctx->tmp[index] = ctx->right[right_index];
            ++right_index;
        }
        ++index;
    }

    while (left_index < ctx->left_size) {
        ctx->tmp[index] = ctx->left[left_index];
        ++left_index;
        ++index;
    }
    while(right_index < ctx->right_size) {
        ctx->tmp[index] = ctx->right[right_index];
        ++right_index;
        ++index;
    }
    if (ctx->tmp != ctx->target) {
        memcpy(ctx->target, ctx->tmp, (ctx->left_size + ctx->right_size) * sizeof(int));
    }
}

void merge_sort(sort_ctx_t* ctx) {
    if (ctx->n <= ctx->m) {
        qsort(ctx->array, ctx->n, sizeof(int), compare);
    } else {
        int mid = ctx->n / 2;
        sort_ctx_t left_ctx = {
            .n = mid,
            .m = ctx->m,
            .array = ctx->array,
            .tmp = ctx->tmp,
        };
        sort_ctx_t right_ctx = {
            .n = ctx->n - mid,
            .m = ctx->m,
            .array = ctx->array + mid,
            .tmp = ctx->tmp + mid,
        };
        merge_sort(&left_ctx);
        merge_sort(&right_ctx);
        merge_ctx_t merge_ctx = {
            .left_size = mid,
            .right_size = ctx->n - mid,
            .left = ctx->array,
            .right = ctx->array + mid,
            .target = ctx->array,
            .tmp = ctx->tmp,
            .m = ctx->m,
        };
        merge(&merge_ctx);
    }
}

int binary_search(int* array, int left, int right, int value) {
    int mid = left + (right - left) / 2;
    if (right - left == 0) {
        return right;
    } else if (array[mid] == value) {
        return mid;
    } else if (array[mid] > value) {
        return binary_search(array, left, mid, value);
    } else {
        return binary_search(array, mid + 1, right, value);
    }
}

void* parallel_merge(void* ptr) {
    merge_ctx_t* ctx = (merge_ctx_t*) ptr;
    if (ctx->threads_num <= 1 || ctx->left_size < ctx->m || ctx->right_size < ctx->m) {
        merge(ctx);
    } else {
        int left_value = ctx->left[ctx->left_size / 2];
        int right_pos = binary_search(ctx->right, 0, ctx->right_size, left_value);

        merge_ctx_t left_ctx = {
            .left_size = ctx->left_size / 2 ,
            .right_size = right_pos,
            .left = ctx->left,
            .right = ctx->right,
            .target = ctx->target,
            .tmp = ctx->tmp,
            .m = ctx->m,
            .threads = ctx->threads + 2,
            .threads_num = (ctx->threads_num - 2) / 2,
        };
        merge_ctx_t right_ctx = {
            .left_size = ctx->left_size - ctx->left_size / 2,
            .right_size = ctx->right_size - right_pos,
            .left = ctx->left + ctx->left_size / 2,
            .right = ctx->right + right_pos,
            .target = ctx->target + ctx->left_size / 2 + right_pos,
            .tmp = ctx->tmp + ctx->left_size / 2 + right_pos,
            .m = ctx->m,
            .threads = ctx->threads + 2 + left_ctx.threads_num,
            .threads_num = ctx->threads_num - 2 - left_ctx.threads_num,
        };
        pthread_t left_thread = *(ctx->threads);
        pthread_t right_thread = *(ctx->threads + 1);
        pthread_create(&left_thread, NULL, parallel_merge, &left_ctx);
        pthread_create(&right_thread, NULL, parallel_merge, &right_ctx);

        pthread_join(left_thread, NULL);
        pthread_join(right_thread, NULL);
    }
    return NULL;
}

void* parallel_merge_sort(void* ptr) {
    sort_ctx_t* ctx = (sort_ctx_t*) ptr;
    if (ctx->threads_num <= 1) {
        merge_sort(ctx);
    } else if (ctx->n <= ctx->m) {
        qsort(ctx->array, ctx->n, sizeof(int), compare);
    } else {
        int mid = ctx->n / 2;
        sort_ctx_t left_ctx = {
            .n = mid,
            .m = ctx->m,
            .array = ctx->array,
            .tmp = ctx->tmp,
            .threads = ctx->threads + 2,
            .threads_num = (ctx->threads_num - 2) / 2,
        };
        sort_ctx_t right_ctx = {
            .n = ctx->n - mid,
            .m = ctx->m,
            .array = ctx->array + mid,
            .tmp = ctx->tmp + mid,
            .threads = ctx->threads + 2 + left_ctx.threads_num,
            .threads_num = ctx->threads_num - 2 - left_ctx.threads_num,
        };
        pthread_t left_thread = *(ctx->threads);
        pthread_t right_thread = *(ctx->threads + 1);
        pthread_create(&left_thread, NULL, parallel_merge_sort, &left_ctx);
        pthread_create(&right_thread, NULL, parallel_merge_sort, &right_ctx);

        pthread_join(left_thread, NULL);
        pthread_join(right_thread, NULL);
        merge_ctx_t merge_ctx = {
            .left_size = mid,
            .right_size = ctx->n - mid,
            .left = ctx->array,
            .right = ctx->array + mid,
            .target = ctx->tmp,
            .tmp = ctx->tmp,
            .m = ctx->m,
        };
        parallel_merge(&merge_ctx);
        memcpy(ctx->array, merge_ctx.target, sizeof(int) * ctx->n);
    }
    return NULL;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        printf("Incorrect number of arguments!");
        return EXIT_FAILURE;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int P = atoi(argv[3]);
    omp_set_num_threads(P);

    int* ar = (int*) malloc(n * sizeof(int));
    int* tmp = (int*) malloc(n * sizeof(int));
    FILE* data = fopen("data.txt", "w");
    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        ar[i] = rand();
        fprintf(data, "%d ", ar[i]);
    }
    fprintf(data, "\n");

    pthread_t threads[P];
    sort_ctx_t ctx = {
        .n = n,
        .m = m,
        .array = ar,
        .tmp = tmp,
        .threads = threads,
        .threads_num = P,
    };

    double start = omp_get_wtime();
    if (P == 1) {
        merge_sort(&ctx);
    } else {
        parallel_merge_sort(&ctx);
    }
    double end = omp_get_wtime();
    double work_time = end - start;
    FILE* stats = fopen("stats.txt", "w");
    fprintf(stats, "%.5fs %d %d %d\n", work_time, n, m, P);
    fclose(stats);

    for (int i = 0; i < n; i++) {
        fprintf(data, "%d ", ctx.array[i]);
    }
    fprintf(data, "\n");
    fclose(data);

    free(ar);
    free(tmp);
    return EXIT_SUCCESS;
}