#include <cstdio>
#include <iostream>
#include "omp.h"

#define N 16
#define K 4
#define NTASKS (N / K)

double A[N][N];
double B[N][N];
bool deps[NTASKS + 1][NTASKS + 1];

void printArray(double tab[N][N])
{
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%4.0lf ", tab[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc, char **argv)
{
  // Initialize the array A[i][j] = i + j in parallel.
  // Using omp for loop
#pragma omp parallel for collapse(2)
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      A[i][j] = i + j;
  printf("Array A:\n");
  printArray(A);
  // TODO: Do the task-based parallel prefix sum
  printf("Array B:\n");
  printArray(B);
  return 0;
}
