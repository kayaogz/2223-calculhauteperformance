#include <iostream>
#include <iomanip>
#include <chrono>
#include "omp.h"

#define P 4

void print(int *m, int N)
{
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
			std::cout << std::setw(5) << m[i*N + j] << " ";
			
		std::cout << std::endl;
	}
}

void multiply(int *m, int N, int *m1, int *m2)
{
	for (int i=0; i<N*N; i++) m[i] = 0;
	
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			for (int k=0; k<N; k++)
				m[i*N + j] += m1[i*N + k] * m2[k*N + j];
}

void add(int *m, int N, int *m1, int *m2, int neg = 1)
{
	for (int i=0; i<N*N; i++)
		m[i] = m1[i] + m2[i] * neg;
}

void partition(int *m, int N, int **ms)
{
	for (int i=0; i<4; i++) ms[i] = new int[N*N/4];
		
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
		{
			int k = (i >= N/2) * 2 + (j >= N/2);
			int ii = i%2;
			int jj = j%2;
			
			ms[k][ii*N/2 + jj] = m[i*N + j];
		}
}

void combine(int *m, int N, int **ms)
{
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
		{
			int k = (i >= N/2) * 2 + (j >= N/2);
			int ii = i%2;
			int jj = j%2;
			
			m[i*N + j] = ms[k][ii*N/2 + jj];
		}
}

void standard_multiply(int *c, int N, int *a, int *b, int K)
{
	if (N <= K) multiply(c, N, a, b);
	else
	{
		int **cs = new int*[4];
		int **as = new int*[4];
		int **bs = new int*[4];
		
		partition(c, N, cs);
		partition(a, N, as);
		partition(b, N, bs);
		
		int **xs = new int*[4];
		int **ys = new int*[4];
		
		for (int i=0; i<4; i++)
		{
			xs[i] = new int[N*N/4];
			ys[i] = new int[N*N/4];
		}
		
		for (int i=0; i<2; i++)
			for (int j=0; j<2; j++)
			{
				#pragma omp task depend(out:xs[i*2 + j])
				standard_multiply(xs[i*2 + j], N/2, as[i*2 + 0], bs[0*2 + j], K);
				
				#pragma omp task depend(out:ys[i*2 + j])
				standard_multiply(ys[i*2 + j], N/2, as[i*2 + 1], bs[1*2 + j], K);
				
				#pragma omp task depend(in:xs[i*2 + j], ys[i*2 + j])
				add(cs[i*2 + j], N/2, xs[i*2 + j], ys[i*2 + j]);
			}
		
		#pragma omp taskwait
		
		combine(c, N, cs);
	}
}

void strassen_multiply(int *c, int N, int *a, int *b, int K)
{
	if (N <= K) multiply(c, N, a, b);
	else
	{
		int **cs = new int*[4];
		int **as = new int*[4];
		int **bs = new int*[4];
		
		partition(c, N, cs);
		partition(a, N, as);
		partition(b, N, bs);
		
		int **ts = new int*[10];
		int **ms = new int*[7];
		int **ns = new int*[4];
		
		for (int i=0; i<10; i++) ts[i] = new int[N*N/4];
		for (int i=0; i<7; i++) ms[i] = new int[N*N/4];
		for (int i=0; i<4; i++) ns[i] = new int[N*N/4];
		
		#pragma omp task depend(out:ts[0])
		add(ts[0], N/2, as[0], as[3]);
		
		#pragma omp task depend(out:ts[1])
		add(ts[1], N/2, bs[0], bs[3]);
		
		#pragma omp task depend(out:ms[0]) depend(in:ts[0], ts[1])
		strassen_multiply(ms[0], N/2, ts[0], ts[1], K);
		
		#pragma omp task depend(out:ts[2])
		add(ts[2], N/2, as[2], as[3]);
		
		#pragma omp task depend(out:ms[1]) depend(in:ts[2])
		strassen_multiply(ms[1], N/2, ts[2], bs[0], K);
		
		#pragma omp task depend(out:ts[3])
		add(ts[3], N/2, bs[1], bs[3], -1);
		
		#pragma omp task depend(out:ms[2]) depend(in:ts[3])
		strassen_multiply(ms[2], N/2, as[0], ts[3], K);
		
		#pragma omp task depend(out:ts[4])
		add(ts[4], N/2, bs[2], bs[0], -1);
		
		#pragma omp task depend(out:ms[3]) depend(in:ts[4])
		strassen_multiply(ms[3], N/2, as[3], ts[4], K);
		
		#pragma omp task depend(out:ts[5])
		add(ts[5], N/2, as[0], as[1]);
		
		#pragma omp task depend(out:ms[4]) depend(in:ts[5])
		strassen_multiply(ms[4], N/2, ts[5], bs[3], K);
		
		#pragma omp task depend(out:ts[6])
		add(ts[6], N/2, as[2], as[0], -1);
		
		#pragma omp task depend(out:ts[7])
		add(ts[7], N/2, bs[0], bs[1]);
		
		#pragma omp task depend(out:ms[5]) depend(in:ts[6], ts[7])
		strassen_multiply(ms[5], N/2, ts[6], ts[7], K);
		
		#pragma omp task depend(out:ts[8])
		add(ts[8], N/2, as[1], as[3], -1);
		
		#pragma omp task depend(out:ts[9])
		add(ts[9], N/2, bs[2], bs[3]);
		
		#pragma omp task depend(out:ms[6]) depend(in:ts[8], ts[9])
		strassen_multiply(ms[6], N/2, ts[8], ts[9], K);
		
		#pragma omp task depend(out:ns[0]) depend(in:ms[0], ms[6])
		add(ns[0], N/2, ms[0], ms[6]);
		
		#pragma omp task depend(out:ns[1]) depend(in:ms[3], ms[4])
		add(ns[1], N/2, ms[3], ms[4], -1);
		
		#pragma omp task depend(in:ns[0], ns[1])
		add(cs[0], N/2, ns[0], ns[1]);
		
		#pragma omp task depend(in:ms[2], ms[4])
		add(cs[1], N/2, ms[2], ms[4]);
		
		#pragma omp task depend(in:ms[1], ms[3])
		add(cs[2], N/2, ms[1], ms[3]);
		
		#pragma omp task depend(out:ns[2]) depend(in:ms[0], ms[1])
		add(ns[2], N/2, ms[0], ms[1], -1);
		
		#pragma omp task depend(out:ns[3]) depend(in:ms[2], ms[5])
		add(ns[3], N/2, ms[2], ms[5]);
		
		#pragma omp task depend(in:ns[2], ns[3])
		add(cs[3], N/2, ns[2], ns[3]);
		
		#pragma omp taskwait
		
		combine(c, N, cs);
	}
}

int main(int argc, char** argv)
{
	if (argc < 3) return 1;
	
	int N = std::atoi(argv[1]);
	int K = std::atoi(argv[2]);
	
	int *a = new int[N*N];
	int *b = new int[N*N];
	int *c = new int[N*N];
	
	for (int i=0; i<N*N; i++) a[i] = b[i] = c[i] = i;
	
	auto t1 = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel num_threads(P)
	#pragma omp single
	standard_multiply(c, N, a, b, K);
	
	auto t2 = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> d1 = t2 - t1;
	std::cout << "Standard method: " << d1.count() << "s" << std::endl;
	
	if (N <= 8)
	{
		std::cout << "Matrix A:" << std::endl;
		print(a, N);
		
		std::cout << "Matrix B:" << std::endl;
		print(b, N);
		
		std::cout << "Matrix C:" << std::endl;
		print(c, N);
		
		std::cout << "--------------------------------------" << std::endl;
	}
	
	auto t3 = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel num_threads(P)
	#pragma omp single
	strassen_multiply(c, N, a, b, K);
	
	auto t4 = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> d2 = t4 - t3;
	std::cout << "Strassen method: " << d2.count() << "s" << std::endl;
	
	if (N <= 8)
	{
		std::cout << "Matrix A:" << std::endl;
		print(a, N);
		
		std::cout << "Matrix B:" << std::endl;
		print(b, N);
		
		std::cout << "Matrix C:" << std::endl;
		print(c, N);
		
		std::cout << "--------------------------------------" << std::endl;
	}
	
	return 0;
}
