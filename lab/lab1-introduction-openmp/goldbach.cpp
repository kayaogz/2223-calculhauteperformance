#include <iostream>
#include <chrono>
#include <cmath>
#include "omp.h"

#define N 8192

/**
  * Calculer si x est un nombre premier ou pas
  * Compute if x is prime number or not
  */
bool estPremier(int x)
{
  if (x < 2) { return false; }
  for (int i = 2; i <= x / 2; i++) {
    if (x % i == 0) { return false; }
  }
  return true;
}

/**
  * Calculer le nombre de paires (i, j) tel que i + j = x et i et j sont des nombres premiers.
  * Compute the number of pairs (i,j) such as i + j = x and i and j are prime numbers.
  */
int goldbach(int x)
{
  int compte = 0;
  if (x <= 2) { return 0; }
  for (int i = 2; i <= x / 2; i++) {
    if (estPremier(i) && estPremier(x - i)) { compte++; }
  }
  return compte;
}

int main()
{
  int goldbachVrai;
  int numPairs[N];
  for (int i = 0; i < N; i++) { numPairs[i] = 0; }

  // La version sequentielle
  // Sequential part
  auto start = std::chrono::high_resolution_clock::now();
  goldbachVrai = 1;
  for (int i = 2; i < N; i += 2) { 
    numPairs[i] = goldbach(i);
    if (numPairs[i] == 0) { goldbachVrai = 0; }
  }
  if (goldbachVrai) { std::cout << "La conjecture de Goldbach est vrai" << std::endl; }
  else { std::cout << "La conjecture de Goldbach est vrai" << std::endl; }
  std::chrono::duration<double> tempsSeq = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps sequentiel: " << tempsSeq.count() << "s\n";

  // Parallelisation sans preciser l'ordonnancement.
  // Parallelise without schedule.
  auto start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  // TODO
  std::chrono::duration<double> temps = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps ordonnancement par defaut: " << temps.count() << "s\n";

  // Parallelisation avec l'ordonnancement static et taille de bloc 256.
  // Parallelise with static schedule and bloc size 256.
  start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  // TODO ...
  temps = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps schedule(static, 256): " << temps.count() << "s\n";

  // Parallelisation avec l'ordonnancement dynamic
  // Parallelise with dynamic schedule.
  start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  // TODO ...
  temps = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps schedule(dynamic): " << temps.count() << "s\n";

  // Parallelisation avec l'ordonnancement dynamic et taille de bloc 256
  // Parallelise with dynamic schedule and bloc size 256.
  start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  // TODO ...
  temps = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps schedule(dynamic, 256): " << temps.count() << "s\n";

  // Parallelisation avec l'ordonnancement guided
  // Parallelise with guided schedule.
  start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  // TODO...
  temps = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps schedule(guided): " << temps.count() << "s\n";

  // Parallelisation avec l'ordonnancement guided et taille de bloc 256
    // Parallelise with guided schedule and bloc size 256.
  start = std::chrono::high_resolution_clock::now();
  // A FAIRE ...
  // TODO ...
  temps = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Temps schedule(guided, 256): " << temps.count() << "s\n";

  return 0;
}
