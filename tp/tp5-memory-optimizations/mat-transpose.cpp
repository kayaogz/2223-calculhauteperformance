#include <iostream>
#include <vector>
#include <chrono>
#include "omp.h"
#include "immintrin.h"

void transAVX8x8_ps(__m256 line[8])
{
  // TODO ...
}

int main()
{
  int N = 1024;
  std::vector<float> A(N * N), B(N * N);
  __m256 tile[8];

  // Initialize the matrix
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) { A[i * N + j] = i - j; }
  }

  // Transpose the matrix sequentially and scalar~(non-vectorized)
  auto start = std::chrono::high_resolution_clock::now();
  // TODO ...
  std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Sequential transpose: " << time.count() << "s\n";

  // Transpose the matrix by 8x8 tiles using AVX transpose
  start = std::chrono::high_resolution_clock::now();
  // TODO ...
  time = std::chrono::high_resolution_clock::now() - start;
  std::cout << "AVX transpose: " << time.count() << "s\n";

  // Transpose the matrix by 8x8 tiles using AVX transpose and in-place transpose
  start = std::chrono::high_resolution_clock::now();
  // TODO ...
  time = std::chrono::high_resolution_clock::now() - start;
  std::cout << "AVX in-place transpose: " << time.count() << "s\n";

  // Transpose the matrix by 8x8 tiles using AVX transpose and in-place transpose and OpenMP tasks
  start = std::chrono::high_resolution_clock::now();
  // TODO ...
  time = std::chrono::high_resolution_clock::now() - start;
  std::cout << "AVX in-place transpose with OpenMP tasks: " << time.count() << "s\n";

  return 0;
}
