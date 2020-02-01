#include "Header.h"




float** matMemAlloc(int numRow, int numCol) //to allocate memomy to a matrix
{
	float** matrix;
	matrix = new float* [numRow];
	for (int i = 0; i < numRow; i++)
	{
		matrix[i] = new float[numCol];
	}
	return matrix;
}



void matMemFree(float** matrix, int numRow) //to free the memory allocated to a matrix
{
	for (int i = 0; i < numRow; i++)
	{
		delete[] matrix[i];
	}
	delete[] matrix;
}



void Initialize(float** A, float** B, int M, int N, int K)
{
	// initializing both matrices using random numbers
	int idx1, idx2;
	float tmp_limit = 0.05f;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		for (idx2 = 0; idx2 < K; idx2++)
		{
			A[idx1][idx2] = -tmp_limit + 2 * tmp_limit * (rand() / float(RAND_MAX));
		}
	}

	for (idx1 = 0; idx1 < K; idx1++)
	{
		for (idx2 = 0; idx2 < N; idx2++)
		{
			B[idx1][idx2] = -tmp_limit + 2 * tmp_limit * (rand() / float(RAND_MAX));
		}
	}
}




void Flatten(float** A, float** B, float* A_flat, float* B_flat_trans, int M, int N, int K)
{
	// flattenning matrix A and flattenning and transposing matrix B
	int idx1, idx2;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		for (idx2 = 0; idx2 < K; idx2++)
		{
			A_flat[idx1 * K + idx2] = A[idx1][idx2];
		}
	}

	for (idx1 = 0; idx1 < N; idx1++)
	{
		for (idx2 = 0; idx2 < K; idx2++)
		{
			B_flat_trans[idx1 * K + idx2] = B[idx2][idx1];
		}
	}
}




void SeqMultiply1(float** A, float** B, float* C, int M, int N, int K)
{
	// sequential matrix multiplication using two-dimensional arrays
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "sequential matrix multiplication using two-dimensional arrays" << endl;

	int idx1, idx2, idx3;
	float tmp_val;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		for (idx2 = 0; idx2 < N; idx2++)
		{
			tmp_val = 0;
			for (idx3 = 0; idx3 < K; idx3++)
			{
				tmp_val += (A[idx1][idx3] * B[idx3][idx2]);
			}
			C[idx1 * N + idx2] = tmp_val;
		}
	}
}




void SeqMultiply2(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	// sequential matrix multiplication, -flatten A, flatten & transposed B
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "sequential matrix multiplication, -flatten A, flatten & transposed B" << endl;

	int idx1, idx2, idx3;
	float tmp_val;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		for (idx2 = 0; idx2 < N; idx2++)
		{
			tmp_val = 0;
			for (idx3 = 0; idx3 < K; idx3++)
			{
				tmp_val += (A_flat[idx1 * K + idx3] * B_flat_trans[idx2 * K + idx3]);
			}
			C[idx1 * N + idx2] = tmp_val;
		}
	}
}




void SeqMultiply3(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* sequential matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "sequential matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                  -pre-calculated indices' offsets" << endl;
	int idx1, idx2, idx3, D1, D2, D3;
	float tmp_val;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			tmp_val = 0; D2 = idx2 * K;
			for (idx3 = 0; idx3 < K; idx3++)
			{
				tmp_val += (A_flat[D1 + idx3] * B_flat_trans[D2 + idx3]);
			}
			C[D3 + idx2] = tmp_val;
		}
	}
}




void SeqMultiply4(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* sequential matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets, -8 times unrolling */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "sequential matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                  -pre-calculated indices' offsets" << endl;
	cout << "                                  -8 times unrolling" << endl;

	int idx1, idx2, idx3, D1, D2, D3;
	float tmp_val, c0, c1, c2, c3, c4, c5, c6, c7;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			tmp_val = 0; D2 = idx2 * K;
			for (idx3 = 0; idx3 < K; idx3+= 8)
			{
				c0 = A_flat[D1 + idx3]     * B_flat_trans[D2 + idx3];
				c1 = A_flat[D1 + idx3 + 1] * B_flat_trans[D2 + idx3 + 1];
				c2 = A_flat[D1 + idx3 + 2] * B_flat_trans[D2 + idx3 + 2];
				c3 = A_flat[D1 + idx3 + 3] * B_flat_trans[D2 + idx3 + 3];
				c4 = A_flat[D1 + idx3 + 4] * B_flat_trans[D2 + idx3 + 4];
				c5 = A_flat[D1 + idx3 + 5] * B_flat_trans[D2 + idx3 + 5];
				c6 = A_flat[D1 + idx3 + 6] * B_flat_trans[D2 + idx3 + 6];
				c7 = A_flat[D1 + idx3 + 7] * B_flat_trans[D2 + idx3 + 7];
				tmp_val += (c0 + c1 + c2 + c3 + c4 + c5 + c6 + c7);
			}
			C[D3 + idx2] = tmp_val;
		}
	}
}




void SeqMultiply5(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* sequential matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets, -SSE intrinsics */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "sequential matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                  -pre-calculated indices' offsets" << endl;
	cout << "                                  -SSE intrinsics" << endl;
	
	int idx1, idx2, idx3, D1, D2, D3;
	__m128 tmp_val, sum;
	int K_iters(K >> 2);
	float result, to_float[4];

	__m128* ptr_A = nullptr;
	__m128* ptr_B = nullptr;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		//std::assume_aligned<ALIGN>(A_flat);
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			sum = _mm_set1_ps(0); D2 = idx2 * K;
			result = 0;
			ptr_A = (__m128*) & A_flat[D1];
			ptr_B = (__m128*) & B_flat_trans[D2];

			for (idx3 = 0; idx3 < K_iters; idx3++, ptr_A++, ptr_B++)
			{
				//tmp_val = _mm_mul_ps(*ptr_A, *ptr_B);
				//sum = _mm_add_ps(tmp_val, sum);
				sum = _mm_fmadd_ps(*ptr_A, *ptr_B, sum);
			}

			_mm_store_ps(to_float, sum);
			for (idx3 = 0; idx3 < 4; idx3++)
			{
				result += to_float[idx3];
			}
			C[D3 + idx2] = result;
		}
	}
}





void SeqMultiply6(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* sequential matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets, -4 times unfolding, -SSE intrinsics */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "sequential matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                  -pre-calculated indices' offsets" << endl;
	cout << "                                  -4 times unfolding" << endl;
	cout << "                                  -SSE intrinsics" << endl;

	int idx1, idx2, idx3, D1, D2, D3;
	__m128 sum, sum0, sum1, sum2, sum3;
	int K_iters(K >> 4);
	float to_float[4];

	__m128* ptr_A = nullptr;
	__m128* ptr_B = nullptr;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			sum = _mm_set1_ps(0); D2 = idx2 * K;
			ptr_A = (__m128*) & A_flat[D1];
			ptr_B = (__m128*) & B_flat_trans[D2];

			for (idx3 = 0; idx3 < K_iters; idx3++)
			{
				sum0 = _mm_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum1 = _mm_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum2 = _mm_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum3 = _mm_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum0 = _mm_add_ps(sum0, sum1);
				sum2 = _mm_add_ps(sum2, sum3);
				sum0 = _mm_add_ps(sum0, sum2);
				sum  = _mm_add_ps(sum0, sum);
			}

			_mm_store_ps(to_float, sum);
			C[D3 + idx2] = to_float[0] + to_float[1] + to_float[2] + to_float[3];
		}
	}
}




void SeqMultiply7(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* sequential matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets, -4 times unfolding, -AVX intrinsics */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "sequential matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                  -pre-calculated indices' offsets" << endl;
	cout << "                                  -4 times unfolding" << endl;
	cout << "                                  -AVX intrinsics" << endl;

	int idx1, idx2, idx3, D1, D2, D3;
	__m256 sum, sum0, sum1, sum2, sum3;
	int K_iters(K >> 5);
	float to_float[8];

	__m256* ptr_A = nullptr;
	__m256* ptr_B = nullptr;

	for (idx1 = 0; idx1 < M; idx1++)
	{
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			sum = _mm256_set1_ps(0); D2 = idx2 * K;
			ptr_A = (__m256*) & A_flat[D1];
			ptr_B = (__m256*) & B_flat_trans[D2];

			for (idx3 = 0; idx3 < K_iters; idx3++)
			{
				sum0 = _mm256_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum1 = _mm256_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum2 = _mm256_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum3 = _mm256_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum0 = _mm256_add_ps(sum0, sum1);
				sum2 = _mm256_add_ps(sum2, sum3);
				sum0 = _mm256_add_ps(sum0, sum2);
				sum  = _mm256_add_ps(sum0, sum);
			}

			_mm256_store_ps(to_float, sum);
			C[D3 + idx2] = to_float[0] + to_float[1] + to_float[2] + 
				to_float[3] + to_float[4] + to_float[5] + to_float[6] + to_float[7];
		}
	}
}




void ParMultiply1(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* parallel matrix multiplication, -flatten A, flatten & transposed B
		-using OpenMP */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "parallel matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                -using OpenMP" << endl;

	int idx1, idx2, idx3;
	float tmp_val;

#pragma omp parallel for private(idx2, idx3, tmp_val)
	for (idx1 = 0; idx1 < M; idx1++)
	{
		for (idx2 = 0; idx2 < N; idx2++)
		{
			tmp_val = 0;
			for (idx3 = 0; idx3 < K; idx3++)
			{
				tmp_val += (A_flat[idx1 * K + idx3] * B_flat_trans[idx2 * K + idx3]);
			}
			C[idx1 * N + idx2] = tmp_val;
		}
	}
}




void ParMultiply2(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* parallel matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets, -using OpenMP */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "parallel matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                -pre-calculated indices' offsets" << endl;
	cout << "                                -using OpenMP" << endl;

	int idx1, idx2, idx3, D1, D2, D3;
	float tmp_val;

#pragma omp parallel for private(idx2, idx3, tmp_val, D1, D2, D3) //schedule(static, 12)
	for (idx1 = 0; idx1 < M; idx1++)
	{
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			tmp_val = 0; D2 = idx2 * K;
			for (idx3 = 0; idx3 < K; idx3++)
			{
				tmp_val += (A_flat[D1 + idx3] * B_flat_trans[D2 + idx3]);
			}
			C[D3 + idx2] = tmp_val;
		}
	}
}




void ParMultiply3(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* parallel matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets, -using OpenMP, -8 times unfolding */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "parallel matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                -pre-calculated indices' offsets" << endl;
	cout << "                                -8 times unfolding" << endl;
	cout << "                                -using OpenMP" << endl;

	int idx1, idx2, idx3, D1, D2, D3;
	float tmp_val, c0, c1, c2, c3, c4, c5, c6, c7;

#pragma omp parallel for private(idx2, idx3, tmp_val, D1, D2, D3, c0, c1, c2, c3, c4, c5, c6, c7) //schedule(static, 12)
	for (idx1 = 0; idx1 < M; idx1++)
	{
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			tmp_val = 0; D2 = idx2 * K;
			for (idx3 = 0; idx3 < K; idx3+=8)
			{
				c0 = A_flat[D1 + idx3]     * B_flat_trans[D2 + idx3];
				c1 = A_flat[D1 + idx3 + 1] * B_flat_trans[D2 + idx3 + 1];
				c2 = A_flat[D1 + idx3 + 2] * B_flat_trans[D2 + idx3 + 2];
				c3 = A_flat[D1 + idx3 + 3] * B_flat_trans[D2 + idx3 + 3];
				c4 = A_flat[D1 + idx3 + 4] * B_flat_trans[D2 + idx3 + 4];
				c5 = A_flat[D1 + idx3 + 5] * B_flat_trans[D2 + idx3 + 5];
				c6 = A_flat[D1 + idx3 + 6] * B_flat_trans[D2 + idx3 + 6];
				c7 = A_flat[D1 + idx3 + 7] * B_flat_trans[D2 + idx3 + 7];
				tmp_val += (c0 + c1 + c2 + c3 + c4 + c5 + c6 + c7);
			}
			C[D3 + idx2] = tmp_val;
		}
	}
}







void ParMultiply4(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* parallel matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets, -using OpenMP, -4 times unfolding 
		-SSE intrinsics */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "parallel matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                -pre-calculated indices' offsets" << endl;
	cout << "                                -4 times unfolding" << endl;
	cout << "                                -SSE intrinsics" << endl;
	cout << "                                -using OpenMP" << endl;

	int idx1, idx2, idx3, D1, D2, D3;
	__m128 sum, sum0, sum1, sum2, sum3;
	int K_iters(K >> 4);
	float result, to_float[4];

	__m128 *ptr_A;
	__m128 *ptr_B;

#pragma omp parallel for private(idx2, idx3, D1, D2, D3, sum, sum0, sum1, sum2, sum3, result, to_float, ptr_A, ptr_B)
	for (idx1 = 0; idx1 < M; idx1++)
	{
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			sum = _mm_set1_ps(0); D2 = idx2 * K;
			result = 0;
			ptr_A = (__m128*) & A_flat[D1];
			ptr_B = (__m128*) & B_flat_trans[D2];

			for (idx3 = 0; idx3 < K_iters; idx3++)
			{
				sum0 = _mm_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum1 = _mm_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum2 = _mm_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum3 = _mm_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum0 = _mm_add_ps(sum0, sum1);
				sum2 = _mm_add_ps(sum2, sum3);
				sum0 = _mm_add_ps(sum0, sum2);
				sum = _mm_add_ps(sum0, sum);
			}

			_mm_store_ps(to_float, sum);
			for (idx3 = 0; idx3 < 4; idx3++)
			{
				result += to_float[idx3];
			}
			C[D3 + idx2] = result;
		}
	}
}





void ParMultiply5(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K)
{
	/* parallel matrix multiplication, -flatten A, flatten & transposed B
		-pre-calculated indices' offsets, -using OpenMP, -4 times unfolding
		-AVX intrinsics */
	cout << "////////////////////////////////////////////////////////////////////" << endl;
	cout << "parallel matrix multiplication, -flatten A, flatten & transposed B" << endl;
	cout << "                                -pre-calculated indices' offsets" << endl;
	cout << "                                -4 times unfolding" << endl;
	cout << "                                -AVX intrinsics" << endl;
	cout << "                                -using OpenMP" << endl;

	int idx1, idx2, idx3, D1, D2, D3;
	__m256 sum, sum0, sum1, sum2, sum3;
	int K_iters(K >> 5);
	float to_float[8];

	__m256* ptr_A;
	__m256* ptr_B;

#pragma omp parallel for private(idx2, idx3, D1, D2, D3, sum, sum0, sum1, sum2, sum3, to_float, ptr_A, ptr_B)
	for (idx1 = 0; idx1 < M; idx1++)
	{
		D1 = idx1 * K; D3 = idx1 * N;
		for (idx2 = 0; idx2 < N; idx2++)
		{
			sum = _mm256_set1_ps(0); D2 = idx2 * K;
			ptr_A = (__m256*) & A_flat[D1];
			ptr_B = (__m256*) & B_flat_trans[D2];

			for (idx3 = 0; idx3 < K_iters; idx3++)
			{
				sum0 = _mm256_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum1 = _mm256_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum2 = _mm256_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum3 = _mm256_mul_ps(*ptr_A, *ptr_B);
				ptr_A++; ptr_B++;

				sum0 = _mm256_add_ps(sum0, sum1);
				sum2 = _mm256_add_ps(sum2, sum3);
				sum0 = _mm256_add_ps(sum0, sum2);
				sum  = _mm256_add_ps(sum0, sum);
			}

			_mm256_store_ps(to_float, sum);
			C[D3 + idx2] = to_float[0] + to_float[1] + to_float[2] + 
				to_float[3] + to_float[4] + to_float[5] + to_float[6] + to_float[7];
		}
	}
}