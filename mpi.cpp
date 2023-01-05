#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <Windows.h>
#include "mpi.h"
#include <clocale>
#include <iostream>
# define K 2000
using namespace std;


// Создание матриц
void create_matrix(double* A, double* B, int N) {
	int i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			A[i * N + j] = B[i * N + j] = i + j;
		}
}

// Инициализация и заполнение матриц
void processInitialize(double*& A, double*& B, double*& C, double*& result, int& N) {

	A = new double[N * N];
	B = new double[N * N];
	C = new double[N * N];
	result = new double[N * N];

	create_matrix(A, B, N);
	for (int i = 0; i < N * N; i++) {
		C[i] = 0;
		result[i] = 0;
	}
}

// умножение матриц
void umn_matr(double* A, double* B, double* C, int N, int nCommRunk, int nCommSize) {
	int i, j, k;
	for (i = 0; i < N; i++) {
		for (j = 0 + nCommRunk; j < N; j += nCommSize) {
			for (k = 0; k < N; k++) {
				C[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}
}


int main(int argc, char* argv[]) {
	setlocale(LC_ALL, "Russian");

	double* A;
	double* B;
	double* C;
	double* result = 0;
	int N = K;
	int nCommRunk, nCommSize, namelen, nCounter;
	int nIntervals;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double t1, t2;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &nCommRunk);
	MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);
	MPI_Get_processor_name(processor_name, &namelen);

	if (nCommRunk == 0) {
		processInitialize(A, B, C, result, N);
		t1 = MPI_Wtime();
		cout << "Начато умножение матриц размерности "<< N <<" элементов с использованием технологии MPI" << endl;
		for (nCounter = 1; nCounter < nCommSize; nCounter++) {
			MPI_Send(&N, 1, MPI_INT, nCounter, 0, MPI_COMM_WORLD);
			MPI_Send(A, N * N, MPI_DOUBLE, nCounter, 1, MPI_COMM_WORLD);
			MPI_Send(B, N * N, MPI_DOUBLE, nCounter, 2, MPI_COMM_WORLD);
			MPI_Send(C, N * N, MPI_DOUBLE, nCounter, 3, MPI_COMM_WORLD);
		}
	}
	else {
		MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		A = new double[N * N];
		B = new double[N * N];
		C = new double[N * N];
		MPI_Recv(A, N * N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(B, N * N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(C, N * N, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	umn_matr(A, B, C, N, nCommRunk, nCommSize);
	MPI_Reduce(C, result, N * N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (nCommRunk == 0) {
		t2 = MPI_Wtime();
		cout << "Время умножения матриц: " << t2 - t1 << endl;
	}
	MPI_Finalize();
}