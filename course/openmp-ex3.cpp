#include <cstdio>
#include "omp.h"

int x;
#pragma omp threadprivate(x)

int main()
{
  int x = 3;
//  printf("%d %x\n", x, &x);
#pragma omp parallel num_threads(2)
  {
    int x;
    x = omp_get_thread_num() + 1000;
    printf("thid = %d, %d %x\n", omp_get_thread_num(), x, &x);
  }
#pragma omp parallel num_threads(2) 
  {
    printf("thid = %d, %d %x\n", omp_get_thread_num(), x, &x);
  }
  return 0;
}
