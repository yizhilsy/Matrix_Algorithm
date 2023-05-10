#pragma once
#include <iostream>
#include "Assistance.h"
using namespace std;
template <typename ElemType>
//矩阵相乘朴素算法
void Basic_Multiplication_numarray(ElemType* a_matrix, ElemType* b_matrix, int a_rows, int a_cols, int b_rows, int b_cols, 
	ElemType* result_matrix)
{
	int result_rows = 0; int result_cols = 0;
	for (int i = 0; i < a_rows; i++)
	{
		for (int j = 0; j < b_cols; j++)
		{
			int sum = 0;
			for (int k = 0; k < a_cols; k++)
			{
				sum = sum + *(a_matrix + i * a_cols + k) * (*(b_matrix + k * b_cols + j));
			}
			*(result_matrix + result_rows * b_cols + j) = sum; result_cols++;
		}
		result_cols = 0; result_rows++;
	}
}
//矩阵相加算法
template <class T>
void Matrix_Add(T* a_matrix, T* b_matrix, int a_rows, int a_cols, int b_rows, int b_cols, T* result_matrix)
{
	if (a_rows != b_rows || a_cols != b_cols) { cout << "Range Error!" << endl; return; }
	else
	{
		for (int i = 0; i < a_rows; i++)
		{
			for (int j = 0; j < a_cols; j++)
			{
				*(result_matrix + i * a_cols + j) = *(a_matrix + i * a_cols + j) + *(b_matrix + i * b_cols + j);
			}
		}
	}
}
//矩阵相减算法
template <class T>
void Matrix_Sub(T* a_matrix, T* b_matrix, int a_rows, int a_cols, int b_rows, int b_cols, T* result_matrix)
{
	if (a_rows != b_rows || a_cols != b_cols) { cout << "Range Error!" << endl; return; }
	else
	{
		for (int i = 0; i < a_rows; i++)
		{
			for (int j = 0; j < b_rows; j++)
			{
				*(result_matrix + i * a_cols + j) = *(a_matrix + i * a_cols + j) - *(b_matrix + i * b_cols + j);
			}
		}
	}
}
template <class T>
void Strassen_Multiplication_numarray(T* a_matrix, T* b_matrix, int a_rows, int a_cols, int b_rows, int b_cols,
	T* result_matrix)
{
	int asize = a_rows > a_cols ? a_rows : a_cols; int bsize = b_rows > b_cols ? b_rows : b_cols;
	int size = asize > bsize ? asize : bsize;
	if (size % 2 == 1) { size++; }
	int* amcopy = new int[size*size]; int* bmcopy = new int[size * size];//开辟堆空间size*size大小
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i < a_rows && j < a_cols) { *(amcopy + i * size + j) = *(a_matrix + i * a_cols + j); }
			else { *(amcopy + i * size + j) = 0; }
			if (i < b_rows && j < b_cols) { *(bmcopy + i * size + j) = *(b_matrix + i * b_cols + j); }
			else { *(bmcopy + i * size + j) = 0; }
		}
	}

	int divide = size / 2;
	int* A11 = new int[divide * divide]; int* A12 = new int[divide * divide]; int* A21 = new int[divide * divide]; int* A22 = new int[divide * divide];
	int* B11 = new int[divide * divide]; int* B12 = new int[divide * divide]; int* B21 = new int[divide * divide]; int* B22 = new int[divide * divide];
	int* C11 = new int[divide * divide]; int* C12 = new int[divide * divide]; int* C21 = new int[divide * divide]; int* C22 = new int[divide * divide];
	int* P1 = new int[divide * divide]; int* P2 = new int[divide * divide]; int* P3 = new int[divide * divide]; int* P4 = new int[divide * divide]; int* P5 = new int[divide * divide];
	int* P6 = new int[divide * divide]; int* P7 = new int[divide * divide];
	int* S1 = new int[divide * divide]; int* S2 = new int[divide * divide]; int* S3 = new int[divide * divide]; int* S4 = new int[divide * divide]; int* S5 = new int[divide * divide];
	int* S6 = new int[divide * divide]; int* S7 = new int[divide * divide]; int* S8 = new int[divide * divide]; int* S9 = new int[divide * divide]; int* S10 = new int[divide * divide];
	for (int i = 0; i < divide; i++)
	{
		for (int j = 0; j < divide; j++)
		{
			*(A11 + i * divide + j) = *(amcopy + i * size + j); *(A12 + i * divide + j) = *(amcopy + i * size + j + divide);
			*(A21 + i * divide + j) = *(amcopy + (i + divide) * size + j); *(A22 + i * divide + j) = *(amcopy + (i + divide) * size + j + divide);
			*(B11 + i * divide + j) = *(bmcopy + i * size + j); *(B12 + i * divide + j) = *(bmcopy + i * size + j + divide);
			*(B21 + i * divide + j) = *(bmcopy + (i + divide) * size + j); *(B22 + i * divide + j) = *(bmcopy + (i + divide) * size + j + divide);
		}
	}
	Matrix_Sub(B12, B22, divide, divide, divide, divide, S1);
	Matrix_Add(A11, A12, divide, divide, divide, divide, S2);
	Matrix_Add(A21, A22, divide, divide, divide, divide, S3);
	Matrix_Sub(B21, B11, divide, divide, divide, divide, S4);
	Matrix_Add(A11, A22, divide, divide, divide, divide, S5);
	Matrix_Add(B11, B22, divide, divide, divide, divide, S6);
	Matrix_Sub(A12, A22, divide, divide, divide, divide, S7);
	Matrix_Add(B21, B22, divide, divide, divide, divide, S8);
	Matrix_Sub(A11, A21, divide, divide, divide, divide, S9);
	Matrix_Add(B11, B12, divide, divide, divide, divide, S10);
	//由于Basic_Multiplication的引用参数，故特设
	Basic_Multiplication_numarray(A11, S1, divide, divide, divide, divide, P1);
	Basic_Multiplication_numarray(S2, B22, divide, divide, divide, divide, P2);
	Basic_Multiplication_numarray(S3, B11, divide, divide, divide, divide, P3);
	Basic_Multiplication_numarray(A22, S4, divide, divide, divide, divide, P4);
	Basic_Multiplication_numarray(S5, S6, divide, divide, divide, divide, P5);
	Basic_Multiplication_numarray(S7, S8, divide, divide, divide, divide, P6);
	Basic_Multiplication_numarray(S9, S10, divide, divide, divide, divide, P7);
	//Calculate matrix:C11
	Matrix_Add(P5, P4, divide, divide, divide, divide, C11); Matrix_Sub(C11, P2, divide, divide, divide, divide, C11);
	Matrix_Add(C11, P6, divide, divide, divide, divide, C11);
	//Calculate matrix:C12
	Matrix_Add(P1, P2, divide, divide, divide, divide, C12);
	//Calculate matrix:C21
	Matrix_Add(P3, P4, divide, divide, divide, divide, C21);
	//Calculate matrix:C22
	Matrix_Add(P5, P1, divide, divide, divide, divide, C22); Matrix_Sub(C22, P3, divide, divide, divide, divide, C22);
	Matrix_Sub(C22, P7, divide, divide, divide, divide, C22);
	for (int i = 0; i < a_rows; i++)
	{
		for (int j = 0; j < b_cols; j++)
		{
			if (i < divide && j < divide){*(result_matrix + i * b_cols + j) = *(C11 + i * divide + j);}
			else if (i < divide && j >= divide){*(result_matrix + i * b_cols + j) = *(C12 + i * divide + j-divide);}
			else if (i >= divide && j < divide){*(result_matrix + i * b_cols + j) = *(C21 + (i-divide) * divide + j);}
			else if(i>=divide&&j>=divide){ *(result_matrix + i * b_cols + j) = *(C22 + (i-divide) * divide + j-divide); }
		}
	}
	//释放堆区数据
	delete[] amcopy; delete[] bmcopy; delete[] A11; delete[] A12; delete[] A21; delete[] A22;
	delete[] B11; delete[] B12; delete[] B21; delete[] B22; delete[] C11; delete[] C12; delete[] C21; delete[] C22;
	delete[] P1; delete[] P2; delete[] P3; delete[] P4; delete[] P5; delete[] P6; delete[] P7;
	delete[] S1; delete[] S2; delete[] S3; delete[] S4; delete[] S5; delete[] S6; delete[] S7; delete[] S8; delete[] S9; delete[] S10;
}
template<class T>
double numarray_determinant_value(T* deterarray,int n)//求行列式函数，运用递归思想
{
	double detval = 0;
	if (n == 1) return *(deterarray);
	T* p_sondet = new T[(n - 1) * (n - 1)];
	for (int i = 0; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			for (int k = 0; k < n ; k++)
			{
				if (k<i)
				{
					*(p_sondet + (j - 1) * (n - 1) + k) = *(deterarray + j * n + k);
				}
				else if (k == i) { continue; }
				else
				{
					*(p_sondet + (j - 1) * (n - 1) + k-1) = *(deterarray + j * n + k);
				}
			}
		}
		detval = detval + *(deterarray + i) * pow(-1.0, i)* numarray_determinant_value(p_sondet,n-1);
	}
	delete[] p_sondet;
	return detval;
}