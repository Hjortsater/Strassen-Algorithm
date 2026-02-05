#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix_alg.h"

#ifndef OPTMULT
#define OPTMULT 1
#endif

#ifndef PRINT
#define PRINT 0
#endif



int main(int argc, char* argv[]){
    if (argc != 2){
        printf("Usage: %s N\n", argv[0]);
        return 1;
    }

    char *end;
    long tmp = strtol(argv[1], &end, 10);
    if (*end != '\0' || tmp <= 0 || tmp > 9999){
        printf("Matrix dimension must be a positive integer between 1 and 9999.\n");
        return 1;
    }

    int N = (int)tmp;
    srand(1);

    Matrix A = matrix_create_rand(N);
    Matrix B = matrix_create_rand(N);
    Matrix C;

    struct timespec wall_start, wall_end;
    clock_gettime(CLOCK_MONOTONIC, &wall_start);

    clock_t cpu_start = clock();

#if OPTMULT
    C = matrix_mul_opt(&A, &B);
#else
    C = matrix_mul_std(&A, &B);
#endif

    clock_t cpu_end = clock();
    clock_gettime(CLOCK_MONOTONIC, &wall_end);

    double wall_sec = (wall_end.tv_sec - wall_start.tv_sec) + 
                      (wall_end.tv_nsec - wall_start.tv_nsec) / 1e9;
    double cpu_sec = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

    printf("\n\n>>>>> wall time: %.6f s\n", wall_sec);
    printf(">>>>> CPU time: %.6f s\n", cpu_sec);



#if PRINT
    printf("Matrix A:\n");
    matrix_print(&A);
    printf("Matrix B:\n");
    matrix_print(&B);
    printf("Matrix C:\n");
    matrix_print(&C);
#endif

    {
        matrix_free(&A);
        matrix_free(&B);
        matrix_free(&C);
    }
    return 0;
}
