#pragma once

#include<stdio.h>
#include<iostream>
#include<vector>
#include<omp.h>
#include<time.h>
#include<emmintrin.h>
#include<immintrin.h>
#include<stdlib.h>
#include<malloc.h>

using namespace std;

#define ALIGN 32
#define K_base 128
#define display(name) printer(#name, (name))


template<typename T> //to print the value and name of a variable
void printer(const char* name, T value) {
	cout << name << " = " << value << endl;
}


float** matMemAlloc(int numRow, int numCol);
void matMemFree(float** matrix, int numRow);

void Initialize(float** A, float** B, int M, int N, int K);

void Flatten(float** A, float** B, float* A_flat, float* B_flat_trans, int M, int N, int K);
void SeqMultiply1(float** A, float** B, float* C, int M, int N, int K);
void SeqMultiply2(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void SeqMultiply3(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void SeqMultiply4(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void SeqMultiply5(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void SeqMultiply6(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void SeqMultiply7(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);

void ParMultiply1(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void ParMultiply2(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void ParMultiply3(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void ParMultiply4(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);
void ParMultiply5(float* A_flat, float* B_flat_trans, float* C, int M, int N, int K);