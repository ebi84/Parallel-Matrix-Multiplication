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

	begin = omp_get_wtime();
	SeqMultiply1(A, B, C, M, N, K);
	end = omp_get_wtime();
	double time1 = (end - begin) * 1000;
	cout << "calc time (ms): " << time1 << endl;

	begin = omp_get_wtime();
	SeqMultiply2(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time2 = (end - begin) * 1000;
	cout << "calc time (ms): " << time2 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time2) / time2 * 100 << endl;

	begin = omp_get_wtime();
	SeqMultiply3(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time3 = (end - begin) * 1000;
	cout << "calc time (ms): " << time3 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time3) / time3 * 100 << endl;

	begin = omp_get_wtime();
	SeqMultiply4(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time4 = (end - begin) * 1000;
	cout << "calc time (ms): " << time4 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time4) / time4 * 100 << endl;

	begin = omp_get_wtime();
	SeqMultiply5(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time5 = (end - begin) * 1000;
	cout << "calc time (ms): " << time5 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time5) / time5 * 100 << endl;

	begin = omp_get_wtime();
	SeqMultiply6(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time6 = (end - begin) * 1000;
	cout << "calc time (ms): " << time6 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time6) / time6 * 100 << endl;

	begin = omp_get_wtime();
	SeqMultiply7(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time7 = (end - begin) * 1000;
	cout << "calc time (ms): " << time7 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time7) / time7 * 100 << endl;

	begin = omp_get_wtime();
	ParMultiply1(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time8 = (end - begin) * 1000;
	cout << "calc time (ms): " << time8 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time8) / time8 * 100 << endl;

	begin = omp_get_wtime();
	ParMultiply2(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time9 = (end - begin) * 1000;
	cout << "calc time (ms): " << time9 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time9) / time9 * 100 << endl;

	begin = omp_get_wtime();
	ParMultiply3(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time10 = (end - begin) * 1000;
	cout << "calc time (ms): " << time10 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time10) / time10 * 100 << endl;

	begin = omp_get_wtime();
	ParMultiply4(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time11 = (end - begin) * 1000;
	cout << "calc time (ms): " << time11 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time11) / time11 * 100 << endl;

	begin = omp_get_wtime();
	ParMultiply5(A_flat, B_flat_trans, C, M, N, K);
	end = omp_get_wtime();
	double time12 = (end - begin) * 1000;
	cout << "calc time (ms): " << time12 << endl;
	cout << "speedup compared to the naive approach(%): " << (time1 - time12) / time12 * 100 << endl;

	_mm_free(C);
	_mm_free(A_flat);
	_mm_free(B_flat_trans);
	matMemFree(A, M);
	matMemFree(B, K);

	return 0;
}