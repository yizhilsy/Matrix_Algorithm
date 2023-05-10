#include <iostream>
#include "Assistance.h"
#include "TriSparseMatrix.h"
#include "Menu.h"
#include "Matrix_Function.h"
using namespace std;
int main()
{
	cout<<"�������������������԰� Version1.0        Arthur:Shiyu Lu"<<endl ;
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
			int n; cout << "��������ʽ����:" << endl; cin >> n;
			cout << "�빹��" <<"n��" << "��A����ʽ:" << endl;
			double* A_matrix = new double[n*n];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++){cin >> *(A_matrix + i * (n)+j);}
			}
			cout<<"����ʽ��ֵΪ:" << numarray_determinant_value(A_matrix, n)<<endl;
			break;
		}
		case 6://��ɢ��ϵ������㵥Ԫ
		{
			cout<<"��ҳ->��ɢ��ѧ���" << endl;
			int a_rows; int a_cols;
			cout << "����Ҫ������A���������:" << endl; cin >> a_rows; cout << "����Ҫ������A���������:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "�빹��" << a_rows << "rows," << a_cols << "cols" << "��A����:" << endl;
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
				cout << "��ҳ->��ɢ��ѧ���->����ϵ���󴫵ݱհ�" << endl;
				Warshall(testtsm01, testresult); testresult.Print();
				break;
			}
			case 2:
			{
				cout << "��ҳ->��ɢ��ѧ���->����ϵ����ԳƱհ�" << endl;
				Relation_Symmetry(testtsm01, testresult); testresult.Print();
				break;
			}
			case 1:
			{
				cout << "��ҳ->��ɢ��ѧ���->����ϵ�����Է��հ�" << endl;
				Relation_Reflect(testtsm01, testresult); testresult.Print();
				break;
			}
			}
			delete[] A_matrix; break;
		}
		case 5:
		{
			cout<<"��ҳ->��������" << endl;
			int n; cout << "����������:" << endl; cin >> n;
			cout << "�빹��" << "n��" << "��A����:" << endl;
			double* A_matrix = new double[n * n];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++) { cin >> *(A_matrix + i * (n)+j); }
			}
			TriSparseMatrix<double> testtsm01(A_matrix, n, n); TriSparseMatrix<double> testresult;
			if (Matrix_Inversion(testtsm01, testresult) == SUCCESS) testresult.Print();
			else cout<<"�þ����޷��������" << endl;
			break;
		}
		case 4://���������˵�Ԫ
		{
			cout << "��ҳ->�������" << endl;
			int a_rows; int a_cols; int b_rows; int b_cols;
			cout << "����Ҫ������A���������:" << endl; cin >> a_rows; cout << "����Ҫ������A���������:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "�빹��" << a_rows << "rows," << a_cols << "cols" << "��A����:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			cout << "����Ҫ������B���������:" << endl; cin >> b_rows; cout << "����Ҫ������B���������:" << endl; cin >> b_cols;
			int* B_matrix = new int[b_rows * b_cols];
			cout << "�빹��" << b_rows << "rows," << b_cols << "cols" << "��B����:" << endl;
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
				cout << "��ҳ->�������->��Ԫ���㷨" << endl;
				Triple_Multiplication(testtsm01, testtsm02, testresult); testresult.Print(); break;
			}
			case 2:
			{
				cout << "��ҳ->�������->Strassen�㷨" << endl;
				Strassen_Multiplication(testtsm01, testtsm02, testresult); testresult.Print(); break;
			}
			case 1:
			{
				cout << "��ҳ->�������->�����㷨" << endl;
				Basic_Multiplication(testtsm01, testtsm02, testresult); testresult.Print(); break;
			}
			}
			delete[]A_matrix; delete[]B_matrix; break;
		}
		case 3://�������ת�õ�Ԫ
		{
			cout << "��ҳ->����ת��" << endl;
			int a_rows; int a_cols;
			cout << "����Ҫ������A���������:" << endl; cin >> a_rows; cout << "����Ҫ������A���������:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "�빹��" << a_rows << "rows," << a_cols << "cols" << "��A����:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			TriSparseMatrix<int> testtsm01(A_matrix, a_rows, a_cols); TriSparseMatrix<int> testresult;
			testtsm01.FastTranspose(testresult); testresult.Print();
			delete[] A_matrix; break;


		}
		case 2://����������Ԫ
		{
			cout << "��ҳ->�������" << endl;
			int a_rows; int a_cols; int b_rows; int b_cols;
			cout << "����Ҫ������A���������:" << endl; cin >> a_rows; cout << "����Ҫ������A���������:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "�빹��" << a_rows << "rows," << a_cols << "cols" << "��A����:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			cout << "����Ҫ������B���������:" << endl; cin >> b_rows; cout << "����Ҫ������B���������:" << endl; cin >> b_cols;
			int* B_matrix = new int[b_rows * b_cols];
			cout << "�빹��" << b_rows << "rows," << b_cols << "cols" << "��B����:" << endl;
			for (int i = 0; i < b_rows; i++)
			{
				for (int j = 0; j < b_cols; j++) { cin >> *(B_matrix + i * b_cols + j); }
			}
			TriSparseMatrix<int> testtsm01(A_matrix, a_rows, a_cols); TriSparseMatrix<int> testtsm02(B_matrix, b_rows, b_cols);
			TriSparseMatrix<int> testresult;
			Matrix_Sub(testtsm01, testtsm02, testresult); testresult.Print();
			delete[]A_matrix; delete[]B_matrix; break;

		}
		case 1://�������ӵ�Ԫ
		{
			cout<<"��ҳ->�������" << endl;
			int a_rows; int a_cols; int b_rows; int b_cols;
			cout << "����Ҫ������A���������:" << endl; cin >> a_rows; cout << "����Ҫ������A���������:" << endl; cin >> a_cols;
			int* A_matrix = new int[a_rows * a_cols];
			cout << "�빹��" << a_rows << "rows," << a_cols << "cols" << "��A����:" << endl;
			for (int i = 0; i < a_rows; i++)
			{
				for (int j = 0; j < a_cols; j++) { cin >> *(A_matrix + i * a_cols + j); }
			}
			cout << "����Ҫ������B���������:" << endl; cin >> b_rows; cout << "����Ҫ������B���������:" << endl; cin >> b_cols;
			int* B_matrix = new int[b_rows * b_cols];
			cout << "�빹��" << b_rows << "rows," << b_cols << "cols" << "��B����:" << endl;
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