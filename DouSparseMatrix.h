#include <iostream>
#include "Assistance.h"
#include "TriSparseMatrix.h"
#include "Double.h"
using namespace std;
template <typename ElemType>
class DouSparseMatrix
{
protected:
	int maxSize;					// ����Ԫ��������
	int rows, cols, num;			// ϡ����������,����������Ԫ����
	Double<ElemType>* douelems;     //��Ԫ��ָ��
	int* row_offsets;               //ÿ�з�0Ԫ�ص���λ��ָ��
public:
	DouSparseMatrix(int rs = DEFAULT_SIZE, int cs = DEFAULT_SIZE, int size = DEFAULT_SIZE);//Ĭ�Ϲ��캯��
	DouSparseMatrix(TriSparseMatrix<ElemType>& tsm);//������Ԫ�鹹���Ԫ�飬ʹ��ʫ���㷨
	~DouSparseMatrix();
	void Print() const;//������ʽ��ӡ,������ʽʹ���˴Ӷ�Ԫ��չ��Ϊ�����ʫ���㷨
};
template<typename ElemType>
DouSparseMatrix<ElemType>::DouSparseMatrix(int rs, int cs, int size)
{
	rows = rs; cols = cs; maxSize = size; douelems = new Double<ElemType>[maxSize]; num = 0;
	row_offsets = new int[rows + 1]; for (int i = 0; i <= rows; i++) { row_offsets[i] = -1; }
}
template<typename ElemType>
DouSparseMatrix<ElemType>::DouSparseMatrix(TriSparseMatrix<ElemType>& tsm)
{
	rows = tsm.GetRows(); cols = tsm.GetCols(); num = tsm.GetNum();
	maxSize = tsm.maxSize; douelems = new Double<ElemType>[maxSize];
	row_offsets = new int[rows + 1];
	if (num > 0)
	{
		int* r_num = new int[rows];
		for (int i = 0; i < rows; i++)
			r_num[i] = 0;
		for (int i = 0; i < num; i++)
			++r_num[tsm.triElems[i].row];
		for (int i = 0; i < rows; i++)
		{
			if (r_num[i] == 0)
			{
				row_offsets[i] = -1;
			}
		}
		int find_1th;
		for (int i = 0; i < rows; i++)
			if (r_num[i] != 0)
			{
				find_1th = i; break;
			}
		int pos_temp = 0; row_offsets[find_1th] = 0; int len_temp = r_num[find_1th];
		for (int i = find_1th + 1; i < rows; i++)
		{
			if (r_num[i] == 0) { continue; }
			else
			{
				row_offsets[i] = pos_temp + len_temp;
				pos_temp = row_offsets[i]; len_temp = r_num[i];
			}
		}
		row_offsets[rows] = num;
		for (int i = 0; i < num; i++)
		{
			douelems[i].col = tsm.triElems[i].col; douelems[i].value = tsm.triElems[i].value;
		}
		delete[] r_num;
	}
}
template<typename ElemType>
void DouSparseMatrix<ElemType>::Print() const
{
	if (num > 0)
	{
		cout << "row_offsets:" << endl;
		for (int i = 0; i < rows + 1; i++)
		{
			cout << row_offsets[i] << " ";
		}
		cout << endl << "��Ԫ�������ʽ:" << endl;
		for (int i = 0; i < num; i++)
		{
			cout << "(" << douelems[i].col << "," << douelems[i].value << ")" << " ";
		}
		cout << endl << "���������ʽ:" << endl;
		int* r_num = new int[rows];
		for (int i = 0; i < rows; i++)
		{
			if (row_offsets[i] == -1)
			{
				r_num[i] = 0;
			}
		}
		int index_1th; int index_2th;
		for (int i = 0; i < rows; i++)
		{
			if (row_offsets[i] != -1)
			{
				index_1th = i; break;
			}
		}
		for (int i = index_1th + 1; i < rows + 1; i++)
		{
			if (row_offsets[i] != -1)
			{
				index_2th = i; break;
			}
		}
		int index_temp1 = row_offsets[index_1th]; int index_temp2 = row_offsets[index_2th];
		for (int i = index_1th; i < rows; i++)
		{
			if (row_offsets[i] == -1)
			{
				continue;
			}
			else
			{
				r_num[i] = index_temp2 - index_temp1;
				if (index_temp2 == num) { break; }
				index_temp1 = index_temp2;
				for (int j = index_2th + 1; j < rows + 1; j++)
				{
					if (row_offsets[j] != -1)
					{
						index_2th = j; index_temp2 = row_offsets[j]; break;
					}
				}
			}
		}
		for (int i = 0; i < rows; i++)
		{
			if (r_num[i] == 0)
			{
				for (int j = 0; j < cols; j++)
				{
					cout << 0 << "\t";
				}
			}
			else
			{
				int endtail; int j = 0;
				for (int ct = 0; ct < r_num[i]; ct++)
				{
					endtail = douelems[row_offsets[i] + ct].col;
					while (j < endtail)
					{
						cout << 0 << "\t"; j++;
					}
					cout << douelems[row_offsets[i] + ct].value << "\t"; j++;
				}
				for (; j < cols; j++) cout << 0 << "\t";
			}
			cout << endl;
		}
		delete[] r_num;
	}
}

template<typename ElemType>
DouSparseMatrix<ElemType>::~DouSparseMatrix()
{
	delete[] douelems; delete[] row_offsets;
}
