#include <chrono>
#include <iostream>
#include <vector>
#include "omp.h"

int main()
{
  int i;
  int N = 10000000;
  std::vector<double> A(N);
  double sum = 0.0;

  auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel
  {
#pragma omp for
    for (i = 0; i < N; i++) {
      A[i] = (double)i;
    } // implicit barrier for all threads
    
    // sum is a shared variable here; reduction clause takes care of race conditions
//#pragma omp for reduction(+:sum)
//    for (i = 0; i < N; i++) {
//      sum += A[i];
//    }

#pragma omp sections
    {
#pragma omp section
      {
        double sumLocale = 0;
        for (int i = 0; i < N / 4; i++) { sumLocale += A[i]; }
#pragma omp atomic
//#pragma omp critical
        sum += sumLocale;
      }
#pragma omp section
      {
        double sumLocale = 0;
        for (int i = N / 4; i < (2 * N) / 4; i++) { sumLocale += A[i]; }
#pragma omp atomic
//#pragma omp critical
        sum += sumLocale;
      }
#pragma omp section
      {
        double sumLocale = 0;
        for (int i = (2 * N) / 4; i < (3 * N) / 4; i++) { sumLocale += A[i]; }
#pragma omp atomic
//#pragma omp critical
        sum += sumLocale;
      }
#pragma omp section
      {
        double sumLocale = 0;
        for (int i = (3 * N) / 4; i < N; i++) { sumLocale += A[i]; }
#pragma omp atomic
//#pragma omp critical
        sum += sumLocale;
      }
    }
  }
  std::cout << "The sum is " << sum << std::endl;
  std::chrono::duration<double> temps = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Time: " << temps.count() << "s\n";

  return 0;
}
