#include "Header.h"



int main()
{
	srand(0);
	int M(4000), N(1200), K(2048);
	if (K % K_base) { cout << "K must be divisible by " << K_base << endl; return -1; }

	float** A = matMemAlloc(M, K);
	float** B = matMemAlloc(K, N);
	float* A_flat       = (float*)_mm_malloc(M * K * sizeof(float), ALIGN);
	float* B_flat_trans = (float*)_mm_malloc(K * N * sizeof(float), ALIGN);
	float* C            = (float*)_mm_malloc(M * N * sizeof(float), ALIGN);

	Initialize(A, B, M, N, K);   // initialize both matrices with random numbers
	Flatten(A, B, A_flat, B_flat_trans, M, N, K);   // flatten matrices

	double begin, end;

	//begin = omp_get_wtime();
	//SeqMultiply1(A, B, C, M, N, K);
	//end = omp_get_wtime();
	//cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	SeqMultiply2(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	SeqMultiply3(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	SeqMultiply4(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	SeqMultiply5(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	SeqMultiply6(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	SeqMultiply7(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	ParMultiply1(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	ParMultiply2(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	ParMultiply3(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	ParMultiply4(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	begin = omp_get_wtime();
	ParMultiply5(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	cout << "calc time (ms): " << (end - begin) * 1000 << endl << endl;

	_mm_free(C);
	_mm_free(A_flat);
	_mm_free(B_flat_trans);
	matMemFree(A, M);
	matMemFree(B, K);

	return 0;
}