#include <iostream>
#include "Assistance.h"
#include "TriSparseMatrix.h"
#include "Menu.h"
#include "Matrix_Function.h"
using namespace std;
int main()
{
	cout<<"轻量化矩阵运算器测试版 Version1.0        Arthur:Shiyu Lu"<<endl ;
	int breakflag = 0;
	while (true)
	{
		menu(); int key1;
		while (true)
		{
			cin >> key1;
			if (key1 >= 0 && key1 <= 7){break;}
			else { cout << "Please Enter the right command number again!" << endl; }
		}
		switch (key1)
		{
		case 7:
		{
			int n; cout << "输入行列式阶数:" << endl; cin >> n;
			cout << "请构建" <<"n阶" << "的A行列式:" << endl;
			double* A_matrix = new double[n*n];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++){cin >> *(A_matrix + i * (n)+j);}
			}
			cout<<"行列式的值为:" << numarray_determinant_value(A_matrix, n)<<endl;
			break;
		}
		case 6://离散关系矩阵计算单元
		{
			cout<<"主页->离散数学板块" << endl;
			int a_rows; int a_cols;
			cout << "输入要构建的A矩阵的行数:" << endl; cin >> a_rows; cout << "输入要构建的A矩阵的列数:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "请构建" << a_rows << "rows," << a_cols << "cols" << "的A矩阵:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			TriSparseMatrix<int> testtsm01(A_matrix, a_rows, a_cols); TriSparseMatrix<int> testresult;
			discrete_math_menu();
			int key3;
			while (true)
			{
				cin >> key3;
				if (key3 >= 1 && key3 <= 3) break;
				else cout << "Please Enter the right command number agagin!" << endl;
			}
			switch (key3)
			{
			case 3:
			{
				cout << "主页->离散数学板块->求解关系矩阵传递闭包" << endl;
				Warshall(testtsm01, testresult); testresult.Print();
				break;
			}
			case 2:
			{
				cout << "主页->离散数学板块->求解关系矩阵对称闭包" << endl;
				Relation_Symmetry(testtsm01, testresult); testresult.Print();
				break;
			}
			case 1:
			{
				cout << "主页->离散数学板块->求解关系矩阵自反闭包" << endl;
				Relation_Reflect(testtsm01, testresult); testresult.Print();
				break;
			}
			}
			delete[] A_matrix; break;
		}
		case 5:
		{
			cout<<"主页->矩阵求逆" << endl;
			int n; cout << "输入矩阵阶数:" << endl; cin >> n;
			cout << "请构建" << "n阶" << "的A矩阵:" << endl;
			double* A_matrix = new double[n * n];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++) { cin >> *(A_matrix + i * (n)+j); }
			}
			TriSparseMatrix<double> testtsm01(A_matrix, n, n); TriSparseMatrix<double> testresult;
			if (Matrix_Inversion(testtsm01, testresult) == SUCCESS) testresult.Print();
			else cout<<"该矩阵无法求逆矩阵！" << endl;
			break;
		}
		case 4://计算矩阵相乘单元
		{
			cout << "主页->矩阵相乘" << endl;
			int a_rows; int a_cols; int b_rows; int b_cols;
			cout << "输入要构建的A矩阵的行数:" << endl; cin >> a_rows; cout << "输入要构建的A矩阵的列数:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "请构建" << a_rows << "rows," << a_cols << "cols" << "的A矩阵:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			cout << "输入要构建的B矩阵的行数:" << endl; cin >> b_rows; cout << "输入要构建的B矩阵的列数:" << endl; cin >> b_cols;
			int* B_matrix = new int[b_rows * b_cols];
			cout << "请构建" << b_rows << "rows," << b_cols << "cols" << "的B矩阵:" << endl;
			for (int i = 0; i < b_rows; i++)
			{
				for (int j = 0; j < b_cols; j++) { cin >> *(B_matrix + i * b_cols + j); }
			}
			TriSparseMatrix<int> testtsm01(A_matrix, a_rows, a_cols); TriSparseMatrix<int> testtsm02(B_matrix, b_rows, b_cols);
			TriSparseMatrix<int> testresult;
			Matrix_Multiplication_menu(); int key2;
			while (true)
			{
				cin >> key2;
				if (key2 >= 1 && key2 <= 3) break;
				else cout << "Please Enter the right command number again!" << endl;
			}
			switch (key2)
			{
			case 3:
			{
				cout << "主页->矩阵相乘->三元组算法" << endl;
				Triple_Multiplication(testtsm01, testtsm02, testresult); testresult.Print(); break;
			}
			case 2:
			{
				cout << "主页->矩阵相乘->Strassen算法" << endl;
				Strassen_Multiplication(testtsm01, testtsm02, testresult); testresult.Print(); break;
			}
			case 1:
			{
				cout << "主页->矩阵相乘->朴素算法" << endl;
				Basic_Multiplication(testtsm01, testtsm02, testresult); testresult.Print(); break;
			}
			}
			delete[]A_matrix; delete[]B_matrix; break;
		}
		case 3://计算矩阵转置单元
		{
			cout << "主页->矩阵转置" << endl;
			int a_rows; int a_cols;
			cout << "输入要构建的A矩阵的行数:" << endl; cin >> a_rows; cout << "输入要构建的A矩阵的列数:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "请构建" << a_rows << "rows," << a_cols << "cols" << "的A矩阵:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			TriSparseMatrix<int> testtsm01(A_matrix, a_rows, a_cols); TriSparseMatrix<int> testresult;
			testtsm01.FastTranspose(testresult); testresult.Print();
			delete[] A_matrix; break;


		}
		case 2://计算矩阵减单元
		{
			cout << "主页->矩阵相减" << endl;
			int a_rows; int a_cols; int b_rows; int b_cols;
			cout << "输入要构建的A矩阵的行数:" << endl; cin >> a_rows; cout << "输入要构建的A矩阵的列数:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "请构建" << a_rows << "rows," << a_cols << "cols" << "的A矩阵:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			cout << "输入要构建的B矩阵的行数:" << endl; cin >> b_rows; cout << "输入要构建的B矩阵的列数:" << endl; cin >> b_cols;
			int* B_matrix = new int[b_rows * b_cols];
			cout << "请构建" << b_rows << "rows," << b_cols << "cols" << "的B矩阵:" << endl;
			for (int i = 0; i < b_rows; i++)
			{
				for (int j = 0; j < b_cols; j++) { cin >> *(B_matrix + i * b_cols + j); }
			}
			TriSparseMatrix<int> testtsm01(A_matrix, a_rows, a_cols); TriSparseMatrix<int> testtsm02(B_matrix, b_rows, b_cols);
			TriSparseMatrix<int> testresult;
			Matrix_Sub(testtsm01, testtsm02, testresult); testresult.Print();
			delete[]A_matrix; delete[]B_matrix; break;

		}
		case 1://计算矩阵加单元
		{
			cout<<"主页->矩阵相加" << endl;
			int a_rows; int a_cols; int b_rows; int b_cols;
			cout << "输入要构建的A矩阵的行数:" << endl; cin >> a_rows; cout << "输入要构建的A矩阵的列数:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "请构建" << a_rows << "rows," << a_cols << "cols" << "的A矩阵:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			cout << "输入要构建的B矩阵的行数:" << endl; cin >> b_rows; cout << "输入要构建的B矩阵的列数:" << endl; cin >> b_cols;
			int* B_matrix = new int[b_rows * b_cols];
			cout << "请构建" << b_rows << "rows," << b_cols << "cols" << "的B矩阵:" << endl;
			for (int i = 0; i < b_rows; i++)
			{
					for (int j = 0; j < b_cols; j++) { cin >> *(B_matrix + i * b_cols + j); }
			}
			TriSparseMatrix<int> testtsm01(A_matrix, a_rows, a_cols); TriSparseMatrix<int> testtsm02(B_matrix, b_rows, b_cols);
			TriSparseMatrix<int> testresult;
			Matrix_Add(testtsm01, testtsm02, testresult); testresult.Print();
			delete[]A_matrix; delete[]B_matrix; break;
		}
		case 0:
		{breakflag = 1; break;	}
		}
		if (breakflag == 1){break;}
	}
	return 0;
}