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

    clock_t t = clock();

#if OPTMULT
    C = matrix_mul_opt(&A, &B);
#else
    C = matrix_mul_std(&A, &B);
#endif
    t = clock() - t;

    printf("\n\n >>>>> mul time: %.6f s\n\n\n", (double)t / CLOCKS_PER_SEC);

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
