#include <stdio.h>
#include <unistd.h>
#include "omp.h"
#include <vector>
#include <assert.h>


int main()
{
  std::vector<int> A(20);
  int x = 3;
#pragma omp parallel num_threads(4)
  {
    int thid = omp_get_thread_num();
    int numth = omp_get_num_threads();
    printf("hello from %d/%d\n", thid, numth);
#pragma omp for
    for (int i = 0; i < A.size(); i++) { 
      A[i] = i;
      printf("i = %d, thid = %d\n", i, thid);
    }
//#pragma omp sections
//    {
//#pragma omp section
//      {
//        for (int i = 0; i < A.size() / 2; i++) { A[i] = i; }
//      }
//#pragma omp section
//      {
//        for (int i = A.size() / 2; i < A.size(); i++) { A[i] = i; }
//      }
//    }
  }
  for (int i = 0; i < A.size(); i++) {
    assert(A[i] == i);
  }

  return 0;
}
