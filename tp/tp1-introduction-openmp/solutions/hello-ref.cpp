#include <cstdio>
#include <iostream>
#include "omp.h"

int main()
{
  omp_set_num_threads(5);
#pragma omp parallel num_threads(7)
  {
    printf("Thread %d/%d\n", omp_get_thread_num(), omp_get_num_threads());
#pragma omp single
    printf("Hello from thread %d\n", omp_get_thread_num());
  }
  return 0;
}
