#ifndef STRASSEN_H
#define STRASSEN_H

#define MAT(p, s, i, j) ((p)[(i)*(s) + (j)])

void strassen_rec(const int *A, const int *B, int *C, int N, 
                  int sA, int sB, int sC, int *ws);

int next_power_of_two(int n);

#endif