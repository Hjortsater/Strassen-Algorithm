#ifndef MATRIX_ALG_H
#define MATRIX_ALG_H

typedef struct {
    int dim;
    int *data;
} Matrix;

Matrix matrix_create_rand(int N);
void matrix_free(Matrix *m);
void matrix_print(const Matrix *m);
Matrix matrix_mul_std(const Matrix *A, const Matrix *B);
Matrix matrix_mul_opt(const Matrix *A, const Matrix *B);

#endif
