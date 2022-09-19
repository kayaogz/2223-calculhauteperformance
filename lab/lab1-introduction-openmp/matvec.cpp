#include <iostream>
#include <chrono>
#include <vector>
#include "omp.h"

#define NREPET 1024

int main(int argc, char **argv)
{
  std::cout << "Matrix-vector product using OpenMP\n";
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " [num-rows/columns]\n";
    std::cout << "  Example: " << argv[0] << " 1024\n";
    return 1;
  }
  int dim = std::atoi(argv[1]); 

  std::vector<double> A(dim * dim);
  std::vector<double> x(dim);
  std::vector<double> b(dim);

  // Initialiser A et x tel que A(i, j) = i + j et b(j) = 1.
  // Initialize A and x such that A(i, j) = i + j and b(j) = 1.
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      A[i * dim + j] = i + j;
    }
    x[i] = 1;
  }

  // Calculer b = A * x NREPET fois en sequentiel
  // Compute b = A * x NREPET times sequentially
  auto start = std::chrono::high_resolution_clock::now();
  for (int repet = 0; repet < NREPET; repet++) {
    // A FAIRE ...
    // TO DO ...
  }
  std::chrono::duration<double> timeSeq = std::chrono::high_resolution_clock::now() - start;
  std::cout << std::scientific << "Sequential execution time: " << timeSeq.count() / NREPET << "s" << std::endl;

  // Calculer b = A * x NREPET fois en parallele avec omp for
  // Compute b = A * x NREPET times in parallel with omp for
  start = std::chrono::high_resolution_clock::now();
  for (int repet = 0; repet < NREPET; repet++) {
    // A FAIRE ...
    // TO DO ...
  }
  std::chrono::duration<double> timePar = std::chrono::high_resolution_clock::now() - start;
  std::cout << std::scientific << "Parallel execution time with omp for: " << timePar.count() / NREPET << "s" <<
    std::endl;

  // Calculer b = A * x NREPET fois en parallele avec omp task
  // Compute b = A * x NREPET times in parallel with omp task
  start = std::chrono::high_resolution_clock::now();
  for (int repet = 0; repet < NREPET; repet++) {
    // A FAIRE ...
    // TO DO ...
  }
  std::chrono::duration<double> timeParTasks = std::chrono::high_resolution_clock::now() - start;
  std::cout << std::scientific << "Parallel execution time with omp tasks: " << timeParTasks.count() / NREPET <<
    "s" << std::endl;

  // Calculer et afficher l'acceleration et l'efficacite de la parallelisation avec omp for
  // Compute and print the acceleration and the parallel efficiency with omp for
  double acceleration = 0.0;
  double efficiency = 0.0;
  // A FAIRE ...
  // TODO
  std::cout << "Acceleration: " << acceleration << std::endl;
  std::cout << "Efficiency: " << efficiency << std::endl;

  // Verifier le resultat. b(i) est cense etre (dim - 1) * dim / 2 + i * dim
  // Verify the result. b(i) is supposed to be (dim - 1) * dim / 2 + i * dim
  for (int i = 0; i < dim; i++) {
    double val = (dim - 1) * (double)dim / 2.0 + (double)i * dim;
    if (b[i] != val) {
      std::cout << "incorrect value: b[" << i << "] = " << b[i] << " != " << val << std::endl;
      break;
    }
    if (i == dim - 1) {
      std::cout << "Matrix-vector product computed successfully!\n";
    }
  }

  return 0;
}
