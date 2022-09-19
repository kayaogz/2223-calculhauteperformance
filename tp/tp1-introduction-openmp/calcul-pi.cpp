#include <chrono>
#include <iostream>
#include "omp.h"

inline double f(double x)
{
  return (4 / (1 + x * x));
}

int main()
{
  int i;
  const int N = 100000000;
  double pi = 0.0;

  // Calculer le pi en sequentiel
  // Sequential compute of pi
  auto start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  // TODO ...
  std::cout << "pi = " << pi << std::endl;
  std::chrono::duration<double> tempsSeq = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps sequentiel: " << tempsSeq.count() << "s\n";

  // Calculer le pi avec omp for et reduction
  // Compute pi with omp for and reduction
  start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  // TODO ...
  std::cout << "pi = " << pi << std::endl;
  std::chrono::duration<double> tempsOmpFor = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps parallel omp for: " << tempsOmpFor.count() << "s\n";

  // Calculer le pi avec boucle parallele faite a la main
  // Compute pi with parallel loop "by hand"
  pi = 0.0;
  start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  std::cout << "pi = " << pi << std::endl;
  std::chrono::duration<double> tempsForMain = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps parallel for a la main: " << tempsForMain.count() << "s\n";

  return 0;
}
