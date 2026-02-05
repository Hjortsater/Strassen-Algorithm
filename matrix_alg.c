#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "matrix_alg.h"
#include "strassen.h"

/*
    Macro to simplify 2D indexing: p=pointer, s=stride(N), i=row, j=col.
    Stride CANNOT be found in the matrix struct defined in the header
    file as its only truly required for the strassen algorithm to work
    I did not want to expose it to main. For the same reason addition
    is defined with stride in mind, and as static.
*/

#define MAT(p, s, i, j) ((p)[(i)*(s) + (j)])

static Matrix _matrix_alloc(int N) {
    Matrix m;
    m.dim = N;
    m.data = malloc(N * N * sizeof(int));
    if (!m.data) {
        printf("Memory allocation failed.\n");
        m.dim = 0;
    }
    return m;
}

Matrix matrix_create_rand(int N){
    Matrix m = _matrix_alloc(N);
    if (!m.data) return m;

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) {
            MAT(m.data, N, i, j) = rand() % 10;
        }
    }
    return m;
}

void matrix_print(const Matrix *m){
    for (int i = 0; i < m->dim; i++){
        for (int j = 0; j < m->dim; j++){
            printf("%5d", MAT(m->data, m->dim, i, j));
        }
        printf("\n");
    }
}

void matrix_free(Matrix *m) {
    if (m->data) {
        free(m->data);
    }
    m->data = NULL;
    m->dim = 0;
}

Matrix matrix_mul_std(const Matrix *A, const Matrix *B) {
    if (A->dim != B->dim) {
        printf("Matrix dimensions must match for multiplication.\n");
        Matrix empty = {0, NULL};
        return empty;
    }

    int N = A->dim;
    Matrix C = _matrix_alloc(N);
    if (!C.data) return C;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            long long sum = 0;
            for (int k = 0; k < N; k++) {
                sum += MAT(A->data, N, i, k) * MAT(B->data, N, k, j);
            }
            MAT(C.data, N, i, j) = sum;
        }
    }

    return C;
}

Matrix matrix_mul_opt(const Matrix *A, const Matrix *B) {
    if (A->dim != B->dim) return (Matrix){0, NULL};

    int N = A->dim;
    Matrix C = _matrix_alloc(N);
    if (!C.data) return C;


    /* Padd the matrices to a power of 2 */

    int P = next_power_of_two(N);

    int *Ap = calloc(P * P, sizeof(int));
    int *Bp = calloc(P * P, sizeof(int));
    int *Cp = calloc(P * P, sizeof(int));

    /* Copy original matrices into padded ones */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            MAT(Ap, P, i, j) = MAT(A->data, N, i, j);
            MAT(Bp, P, i, j) = MAT(B->data, N, i, j);
        }

    int *workspace = calloc(32 * P * P, sizeof(int));

    #pragma omp parallel
    {
        #pragma omp single
        {
            strassen_rec(Ap, Bp, Cp, P, P, P, P, workspace);
        }
    }


    /* Unpad result */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            MAT(C.data, N, i, j) = MAT(Cp, P, i, j);

    free(Ap);
    free(Bp);
    free(Cp);
    free(workspace);

    
    return C;
}

