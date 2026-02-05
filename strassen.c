#include "strassen.h"

#ifndef OPTCUTOFF
#define OPTCUTOFF 32
#endif

#ifndef PARALLEL_CUTOFF
#define PARALLEL_CUTOFF 64
#endif


/* Static helpers: only strassen.c can see these */
static void _mat_add(const int *A, const int *B, int *C, int n, int sA, int sB, int sC) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            MAT(C, sC, i, j) = MAT(A, sA, i, j) + MAT(B, sB, i, j);
}

static void _mat_sub(const int *A, const int *B, int *C, int n, int sA, int sB, int sC) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            MAT(C, sC, i, j) = MAT(A, sA, i, j) - MAT(B, sB, i, j);
}

int next_power_of_two(int n) {
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

void strassen_rec(const int *A, const int *B, int *C, int N, int sA, int sB, int sC, int *ws) {
    /* Base Case: For small matrices, use the standard O(n^3) approach */
    if (N <= OPTCUTOFF) {
        /* cache-friendly i-k-j ordering: zero C then accumulate by k */
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                MAT(C, sC, i, j) = 0;

        for (int i = 0; i < N; i++) {
            for (int k = 0; k < N; k++) {
                int a_ik = MAT(A, sA, i, k);
                for (int j = 0; j < N; j++) {
                    MAT(C, sC, i, j) += a_ik * MAT(B, sB, k, j);
                }
            }
        }
        return;
    }

    
    /* 1. Calculate dimensions and prepare workspace pointers */
    int half = (N >> 1) + (N & 1);
    int size = half * half;
    {    
        int *m1 = ws;
        int *m2 = ws + size;
        int *m3 = ws + size * 2;
        int *m4 = ws + size * 3;
        int *m5 = ws + size * 4;
        int *m6 = ws + size * 5;
        int *m7 = ws + size * 6;
        int *t1 = ws + size * 7;
        int *t2 = ws + size * 8;

        int *next_ws = ws + size * 9;

        int ws_stride = size * 9;

        int *ws1 = next_ws + ws_stride * 0;
        int *ws2 = next_ws + ws_stride * 1;
        int *ws3 = next_ws + ws_stride * 2;
        int *ws4 = next_ws + ws_stride * 3;
        int *ws5 = next_ws + ws_stride * 4;
        int *ws6 = next_ws + ws_stride * 5;
        int *ws7 = next_ws + ws_stride * 6;

        /* 2. Subdivide input matrices using pointers and stride */
        const int *A11 = A;
        const int *A12 = A + half;
        const int *A21 = A + half * sA;
        const int *A22 = A + half * sA + half;

        const int *B11 = B;
        const int *B12 = B + half;
        const int *B21 = B + half * sB;
        const int *B22 = B + half * sB + half;

        int *C11 = C;
        int *C12 = C + half;
        int *C21 = C + half * sC;
        int *C22 = C + half * sC + half;

        /* 3. Compute M1 through M7 */
        
        _mat_add(A11, A22, t1, half, sA, sA, half);
        _mat_add(B11, B22, t2, half, sB, sB, half);
        #pragma omp task shared(m1) if (N > PARALLEL_CUTOFF)
        strassen_rec(t1, t2, m1, half, half, half, half, ws1);

        // M2 = (A21 + A22) * B11
        _mat_add(A21, A22, t1, half, sA, sA, half);
        #pragma omp task shared(m2) if (N > PARALLEL_CUTOFF)
        strassen_rec(t1, B11, m2, half, half, sB, half, ws2);

        // M3 = A11 * (B12 - B22)
        _mat_sub(B12, B22, t2, half, sB, sB, half);
        #pragma omp task shared(m3) if (N > PARALLEL_CUTOFF)
        strassen_rec(A11, t2, m3, half, sA, half, half, ws3);

        // M4 = A22 * (B21 - B11)
        _mat_sub(B21, B11, t2, half, sB, sB, half);
        #pragma omp task shared(m4) if (N > PARALLEL_CUTOFF)
        strassen_rec(A22, t2, m4, half, sA, half, half, ws4);

        // M5 = (A11 + A12) * B22
        _mat_add(A11, A12, t1, half, sA, sA, half);
        #pragma omp task shared(m5) if (N > PARALLEL_CUTOFF)
        strassen_rec(t1, B22, m5, half, half, sB, half, ws5);

        // M6 = (A21 - A11) * (B11 + B12)
        _mat_sub(A21, A11, t1, half, sA, sA, half);
        _mat_add(B11, B12, t2, half, sB, sB, half);
        #pragma omp task shared(m6) if (N > PARALLEL_CUTOFF)
        strassen_rec(t1, t2, m6, half, half, half, half, ws6);

        // M7 = (A12 - A22) * (B21 + B22)
        _mat_sub(A12, A22, t1, half, sA, sA, half);
        _mat_add(B21, B22, t2, half, sB, sB, half);
        #pragma omp task shared(m7) if (N > PARALLEL_CUTOFF)
        strassen_rec(t1, t2, m7, half, half, half, half, ws7);

        #pragma omp taskwait

        /* 4. Combine M-matrices into the quadrants of C */

        // C11 = M1 + M4 - M5 + M7
        _mat_add(m1, m4, t1, half, half, half, half);
        _mat_sub(t1, m5, t2, half, half, half, half);
        _mat_add(t2, m7, C11, half, half, half, sC);

        // C12 = M3 + M5
        _mat_add(m3, m5, C12, half, half, half, sC);

        // C21 = M2 + M4
        _mat_add(m2, m4, C21, half, half, half, sC);

        // C22 = M1 - M2 + M3 + M6
        _mat_sub(m1, m2, t1, half, half, half, half);
        _mat_add(t1, m3, t2, half, half, half, half);
        _mat_add(t2, m6, C22, half, half, half, sC);
    }
}
