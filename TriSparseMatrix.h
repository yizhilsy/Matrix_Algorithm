#ifndef __TRI_SPARSE_MATRIX_H__
#define __TRI_SPARSE_MATRIX_H__

#include "Assistance.h"				// ���������
#include "Triple.h"				// ��Ԫ����
#include "Matrix_Function.h"//����������ļ�
template <typename ElemType>	class DouSparseMatrix;

// ϡ�������Ԫ��˳�����
template<class ElemType>
class TriSparseMatrix
{
protected:
	// ϡ�������Ԫ��˳�������ݳ�Ա:
	Triple<ElemType>* triElems;		// �洢ϡ��������Ԫ���
	int maxSize;					// ����Ԫ��������
	int rows, cols, num;			// ϡ����������,����������Ԫ����

public:
	// ϡ�������Ԫ��˳���ĺ�����Ա�� 
	TriSparseMatrix(int rs = DEFAULT_SIZE, int cs = DEFAULT_SIZE, int size = DEFAULT_SIZE);// ����һ��rs��cs�з���Ԫ��������Ϊsize�Ŀ�ϡ�����
	TriSparseMatrix(ElemType* numarray,int rs,int cs , int size = DEFAULT_SIZE);	//�Ӷ�ά���鹹����Ԫ���
	~TriSparseMatrix();				// ��������
	int GetRows() const;			// ����ϡ���������
	int GetCols() const;			// ����ϡ���������
	int GetNum() const;				// ����ϡ��������Ԫ����
	Status SetElem(int r, int c, const ElemType& v);	// ����ָ��λ�õ�Ԫ��ֵ
	Status GetElem(int r, int c, ElemType& v);			// ��ָ��λ�õ�Ԫ��ֵ
	TriSparseMatrix(const TriSparseMatrix<ElemType>& copy);	// ���ƹ��캯��
	TriSparseMatrix<ElemType>& operator =(const TriSparseMatrix<ElemType>& copy);// ��ֵ���������
	void SimpleTranspose(TriSparseMatrix<ElemType>& b);// ϡ�����ļ�ת���㷨
	void FastTranspose(TriSparseMatrix<ElemType>& b);// ϡ�����Ŀ���ת���㷨
	ElemType* Unseal();//���⺯������ѹ������Ԫ���ڶ�����ѹ�󿽱�һ�ݲ����ص�ַ
	void Print() const;//��ӡ����Ԫ��ϡ�����
	friend class DouSparseMatrix<ElemType>;//��Ԫ��������Ԫ�����һ������
	//����/��ɢ��Ԫ��������
	template<class T>
	friend Status Basic_Multiplication(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result);
	template<class T>
	friend Status Strassen_Multiplication(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result);
	template<class T>
	friend Status Triple_Multiplication(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result);
	template<class T>
	friend Status Matrix_Add(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result);
	template<class T>
	friend Status Matrix_Sub(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result);
	template<class T>
	friend Status Relation_Reflect(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& result);
	template<class T>
	friend Status Relation_Symmetry(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& result);
	template<class T>
	friend Status Warshall(TriSparseMatrix<T>& tsma,TriSparseMatrix<T>& result);
	template<class T>
	friend Status Matrix_Inversion(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& result);
};

// ϡ�������Ԫ��˳������ʵ�ֲ���
template <class ElemType>
TriSparseMatrix<ElemType>::TriSparseMatrix(int r, int c, int size)
// ��������� ����һ��r��c�з���Ԫ��������Ϊsize�Ŀ�ϡ�����
{
	if (r < 1 || c < 1)
		throw Error("������������Ч!");	// �׳��쳣
	maxSize = size;						// ����Ԫ��������
	rows = r;							// ����
	cols = c;							// ����
	num = 0;							// ����Ԫ�ظ���
	triElems = new Triple<ElemType>[maxSize];	// ����洢�ռ�
}

template <class ElemType>
TriSparseMatrix<ElemType>::TriSparseMatrix(ElemType* numarray,int rs, int cs, int size)
{
	rows = rs; cols = cs; maxSize = size; this->triElems=  new Triple<ElemType>[maxSize]; num = 0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (*(numarray + i * cols + j) != 0) { this->SetElem(i, j, *(numarray + i * cols + j)); }
		}
	}
}

template <class ElemType>
TriSparseMatrix<ElemType>::~TriSparseMatrix()
// ���������ϡ�������ռ�ÿռ�
{
	if (triElems != NULL) delete[]triElems; // �ͷŴ洢�ռ�
}

template <class ElemType>
int TriSparseMatrix<ElemType>::GetRows() const
// �������������ϡ���������
{
	return rows;					// ��������
}

template <class ElemType>
int TriSparseMatrix<ElemType>::GetCols() const
// �������������ϡ���������
{
	return cols;					// ��������
}

template <class ElemType>
int TriSparseMatrix<ElemType>::GetNum() const
// �������������ϡ��������Ԫ����
{
	return num;						// ���ط���Ԫ����
}

template <class ElemType>
Status TriSparseMatrix<ElemType>::SetElem(int r, int c, const ElemType& v)
// �������������±귶Χ��,�򷵻�RANGE_ERROR,������,�򷵻�OVER_FLOW,����
//	��SUCCESS
{
	if (r >= rows || c >= cols || r < 0 || c < 0)
		return RANGE_ERROR;					// �±귶Χ��

	int i, j;								// ��������
	for (j = num - 1; j >= 0 &&
		(r < triElems[j].row || r == triElems[j].row && c < triElems[j].col); j--);// ������Ԫ��λ��

	if (j >= 0 && triElems[j].row == r && triElems[j].col == c) {	// �ҵ���Ԫ��
		if (v == 0) {	// ɾ����Ԫ��
			for (i = j + 1; i < num; i++)
				triElems[i - 1] = triElems[i];	// ǰ�ƴ�j+1��ʼ����Ԫ��
			num--;								// ɾ����Ԫ���,����Ԫ�����Լ�1
		}
		else	// �޸�Ԫ��ֵ
			triElems[j].value = v;
		return SUCCESS;						// �ɹ�
	}
	else if (v != 0) {
		if (num < maxSize) {	// ����Ԫ��(r, c, v)���뵽��Ԫ�����
			for (i = num - 1; i > j; i--)	// ����Ԫ��	
				triElems[i + 1] = triElems[i];
			// j + 1Ϊ�ճ��Ĳ���λ��
			triElems[j + 1].row = r;		// ��
			triElems[j + 1].col = c;		// ��
			triElems[j + 1].value = v;		// ����Ԫ��ֵ
			num++;							// ������Ԫ���,����Ԫ�����Լ�1
			return SUCCESS;					// �ɹ�
		}
		else	// ���
			return OVER_FLOW;				// ���ʱ����OVER_FLOW
	}
	return SUCCESS;							// �ɹ�
}

template <class ElemType>
Status TriSparseMatrix<ElemType>::GetElem(int r, int c, ElemType& v)
// �������������±귶Χ��,�򷵻�RANGE_ERROR,���򷵻�SUCCESS,����v����ָ��λ��Ԫ��ֵ
{
	if (r >= rows || c >= cols || r < 0 || c < 0)
		return RANGE_ERROR;			// �±귶Χ��

	int j;							// ��������


	for (j = num - 1; j >= 0 &&
		(r < triElems[j].row || r == triElems[j].row && c < triElems[j].col); j--);// ����ָ��λ�õ���Ԫ��

	if (j >= 0 && triElems[j].row == r && triElems[j].col == c)	// �ҵ���Ԫ��
		v = triElems[j].value;		// ��v����ָ��λ��Ԫ��ֵ
	else	// δ�ҵ���Ԫ��
		v = 0;						// δ�ҵ���Ԫ��,��ʾ0Ԫ��ֵ
	return SUCCESS;					// �ɹ�
}

template <class ElemType>
TriSparseMatrix<ElemType>::TriSparseMatrix(const TriSparseMatrix<ElemType>& copy)
// �����������ϡ�����copy������ϡ����󡪡����ƹ��캯��
{
	maxSize = copy.maxSize;							// ������Ԫ�ظ���
	triElems = new Triple<ElemType>[maxSize];		// ����洢�ռ�
	rows = copy.rows;								// ��������
	cols = copy.cols;								// ��������
	num = copy.num;									// ���Ʒ���Ԫ�ظ���
	triElems = new Triple<ElemType>[maxSize];		// Ϊ��Ԫ�����洢�ռ�
	for (int i = 0; i < num; i++)	// ������Ԫ��
		triElems[i] = copy.triElems[i];
}

template <class ElemType>
TriSparseMatrix<ElemType>& TriSparseMatrix<ElemType>::operator =(const TriSparseMatrix<ElemType>& copy)
// �����������ϡ�����copy��ֵ����ǰϡ����󡪡���ֵ���������
{
	if (&copy != this) {
		maxSize = copy.maxSize;						// ������Ԫ�ظ���
		if (triElems != NULL) delete[]triElems;	// �ͷŴ洢�ռ�
		triElems = new Triple<ElemType>[maxSize];	// ����洢�ռ�
		rows = copy.rows;							// ��������
		cols = copy.cols;							// ��������
		num = copy.num;								// ���Ʒ���Ԫ�ظ���

		for (int i = 0; i < num; i++)	// ������Ԫ��
			triElems[i] = copy.triElems[i];
	}
	return *this;
}

template<class ElemType>
void TriSparseMatrix<ElemType>::SimpleTranspose(TriSparseMatrix<ElemType>& b)
// ���������ϡ�����ļ�ת���㷨�����������Ԫ��˳���b�� 
{
	b.rows = cols;
	b.cols = rows;
	b.num = num;
	b.maxSize = maxSize;
	delete[]b.triElems;
	b.triElems = new Triple<ElemType>[b.maxSize];

	if (b.num > 0) {
		int i = 0;						// ϡ�����b����һ����Ԫ��Ĵ��λ��
		for (int col = 0; col < cols; col++)
			for (int j = 0; j < num; j++)	// ����a�����е�col�е���Ԫ��
				if (triElems[j].col == col) {
					b.triElems[i].row = triElems[j].col;
					b.triElems[i].col = triElems[j].row;
					b.triElems[i].value = triElems[j].value;
					i++;
				}
	}
}

template<class ElemType>
void TriSparseMatrix<ElemType>::FastTranspose(TriSparseMatrix<ElemType>& b)//�����㷨��
// ���������ϡ�����Ŀ���ת���㷨�����������Ԫ��˳���b�� 
{
	b.rows = cols;
	b.cols = rows;
	b.num = num;
	b.maxSize = maxSize;
	delete[]b.triElems;
	b.triElems = new Triple<ElemType>[b.maxSize];

	int* cNum = new int[cols + 1];	// ���ԭ������ÿһ�еķ���Ԫ����
	int* cPos = new int[cols + 1];	// ���ԭ������ÿһ�еĵ�һ������Ԫ��b�еĴ洢λ��
	int col;
	int i;

	if (b.num > 0) {
		for (col = 0; col < cols; col++) cNum[col] = 0;	// ��ʼ��cNum
		for (i = 0; i < num; i++)
			++cNum[triElems[i].col];		// ͳ��ԭ������ÿһ�еķ���Ԫ����
		cPos[0] = 0;						// ��һ�еĵ�һ������Ԫ��b�洢����ʼλ��
		for (col = 1; col < cols; col++)	// ѭ����ÿһ�еĵ�һ������Ԫ��b�洢����ʼλ��
			cPos[col] = cPos[col - 1] + cNum[col - 1];

		for (i = 0; i < num; i++) {	// ѭ������ԭ�������Ԫ��
			int j = cPos[triElems[i].col];
			// ���ڱ�ʾb��ǰ�е���һ������Ԫ��Ԫ��Ĵ洢λ��
			b.triElems[j].row = triElems[i].col;
			b.triElems[j].col = triElems[i].row;
			b.triElems[j].value = triElems[i].value;
			++cPos[triElems[i].col];
			// b��ǰ�е���һ������Ԫ��Ԫ��Ĵ洢��λ��
		}
	}
	delete[]cNum;
	delete[]cPos;
}

template<class ElemType>
ElemType* TriSparseMatrix<ElemType>::Unseal()
{
	ElemType* point = new ElemType[rows * cols];
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++){*(point + i * cols + j) = 0;}
	}
	for (int i = 0; i < num; i++)
	{
		*(point + triElems[i].row * cols + triElems[i].col) = triElems[i].value;
	}
	return point;//�����ڶ������ٿռ�ָ�룬�ں�����ʹ�ú���Ҫ�ǵ�ɾ��
}

template<class ElemType>
void TriSparseMatrix<ElemType>::Print() const
{
	cout << "��Ԫ�������ʽ:" << endl;
	for (int i = 0; i < num; i++)
	{
		cout << "(" << triElems[i].row << "," << triElems[i].col << "," << triElems[i].value << ");";
	}
	cout << endl << "���������ʽ:" << endl; int flag = 0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			for (int ct = 0; ct < num; ct++)
			{
				if (triElems[ct].row == i && triElems[ct].col == j)
				{
					cout << triElems[ct].value << "\t"; flag = 1; break;
				}
			}
			if (flag == 0) { cout << 0 << "\t"; }flag = 0;
		}
		cout << endl;
	}
}

//����/��ɢ��Ԫ��������
template <class T>
Status Basic_Multiplication(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result)
{
	if (tsma.cols != tsmb.rows) { return FAIL; }
	else
	{
		T* p1 = tsma.Unseal(); T* p2 = tsmb.Unseal();
		result.rows = tsma.rows; result.cols = tsmb.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
		result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
		T* presult = new T[tsma.rows * tsmb.cols];
		Basic_Multiplication_numarray(p1, p2, tsma.rows, tsma.cols, tsmb.rows, tsmb.cols, presult);
		for (int i = 0; i < tsma.rows; i++)
		{
			for (int j = 0; j < tsmb.cols; j++)
			{
				if (*(presult + i * tsmb.cols + j) != 0)	result.SetElem(i, j, *(presult + i * tsmb.cols + j));
			}
		}
		delete[] p1; delete[] p2; delete[] presult;//ɾ���ڶ�����ѹ��Ԫ�� 
		return SUCCESS;
	}
}

template <class T>
Status Strassen_Multiplication(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result)
{
	if (tsma.cols != tsmb.rows) { return FAIL; }
	else
	{
		T* p1 = tsma.Unseal(); T* p2 = tsmb.Unseal();
		result.rows = tsma.rows; result.cols = tsmb.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
		result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
		T* presult = new T[tsma.rows * tsmb.cols];
		Strassen_Multiplication_numarray(p1, p2, tsma.rows, tsma.cols, tsmb.rows, tsmb.cols, presult);
		for (int i = 0; i < tsma.rows; i++)
		{
			for (int j = 0; j < tsmb.cols; j++)
			{
				if (*(presult + i * tsmb.cols + j) != 0)	result.SetElem(i, j, *(presult + i * tsmb.cols + j));
			}
		}
		delete[] p1; delete[] p2; delete[] presult;//ɾ���ڶ�����ѹ��Ԫ�� 
		return SUCCESS;
	}
}

template<class T>
Status Triple_Multiplication(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result)
{
	if (tsma.cols != tsmb.rows) { return FAIL; }
	else
	{
		T* p2 = tsmb.Unseal();
		result.rows = tsma.rows; result.cols = tsmb.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
		result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
		//����row_offsets�Լ�r_num�����������ĳ˻�
		int* row_offsets = new int[tsma.rows + 1];int* r_num = new int[tsma.rows];
		for (int i = 0; i < tsma.rows; i++)
			r_num[i] = 0;
		for (int i = 0; i < tsma.num; i++)
			++r_num[tsma.triElems[i].row];
		for (int i = 0; i < tsma.rows; i++)
		{
			if (r_num[i] == 0)
			{
				row_offsets[i] = -1;
			}
		}
		int find_1th;
		for (int i = 0; i < tsma.rows; i++)
			if (r_num[i] != 0){find_1th = i; break;}
		int pos_temp = 0; row_offsets[find_1th] = 0; int len_temp = r_num[find_1th];
		for (int i = find_1th + 1; i < tsma.rows; i++)
		{
			if (r_num[i] == 0) { continue; }
			else
			{
				row_offsets[i] = pos_temp + len_temp;
				pos_temp = row_offsets[i]; len_temp = r_num[i];
			}
		}
		row_offsets[tsma.rows] = tsma.num;
		//����row_offsets�Լ�r_num�������ĳ˻�
		int sum = 0;
		for (int i = 0; i < tsma.rows; i++)
		{
			for (int j = 0; j < tsmb.cols; j++)
			{
				sum = 0;
				for (int k = row_offsets[i]; k < row_offsets[i] + r_num[i]; k++)
				{
					sum = sum + tsma.triElems[k].value * *(p2 + tsma.triElems[k].col * tsmb.cols + j); 
				}
				result.SetElem(i, j, sum);
			}
		}
		delete[] p2; //ɾ���ڶ�����ѹ��Ԫ�� 
		delete[] r_num; delete[] row_offsets;
		return SUCCESS;
	}
}

template<class T>
Status Matrix_Add(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result)
{
	if (tsma.rows != tsmb.rows || tsma.cols != tsmb.cols)
	{
		return FAIL;
	}
	else
	{
		T* p1 = tsma.Unseal(); T* p2 = tsmb.Unseal();
		result.rows = tsma.rows; result.cols = tsmb.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
		result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
		T* presult = new T[tsma.rows * tsma.cols];
		Matrix_Add(p1, p2, tsma.rows, tsma.cols, tsmb.rows, tsmb.cols, presult);
		for (int i = 0; i < tsma.rows; i++)
		{
			for (int j = 0; j < tsma.cols; j++)
			{
				if (*(presult + i * tsma.cols + j) != 0)	result.SetElem(i, j, *(presult + i * tsma.cols + j));
			}
		}
		delete[] p1; delete[] p2; delete[] presult;//ɾ���ڶ�����ѹ��Ԫ�� 
		return SUCCESS;
	}
}

template<class T>
Status Matrix_Sub(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& tsmb, TriSparseMatrix<T>& result)
{
	if (tsma.rows != tsmb.rows || tsma.cols != tsmb.cols)
	{
		return FAIL;
	}
	else
	{
		T* p1 = tsma.Unseal(); T* p2 = tsmb.Unseal();
		result.rows = tsma.rows; result.cols = tsmb.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
		result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
		T* presult = new T[tsma.rows * tsma.cols];
		Matrix_Sub(p1, p2, tsma.rows, tsma.cols, tsmb.rows, tsmb.cols, presult);
		for (int i = 0; i < tsma.rows; i++)
		{
			for (int j = 0; j < tsma.cols; j++)
			{
				if (*(presult + i * tsma.cols + j) != 0)	result.SetElem(i, j, *(presult + i * tsma.cols + j));
			}
		}
		delete[] p1; delete[] p2; delete[] presult;//ɾ���ڶ�����ѹ��Ԫ�� 
		return SUCCESS;
	}
}

template<class T>
Status Relation_Reflect(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& result)
{
	result.rows = tsma.rows; result.cols = tsma.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
	result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
	T* IA=new T[tsma.rows*tsma.cols];
	for (int i = 0; i < tsma.rows; i++)
	{
		for (int j = 0; j < tsma.cols; j++)
		{
			if (i == j) *(IA+i*tsma.cols+j) = 1;
			else *(IA + i * tsma.cols + j) = 0;
		}
	}
	TriSparseMatrix<T> IA_Matrix(IA,tsma.rows,tsma.cols);
	Matrix_Add(tsma, IA_Matrix, result); delete[] IA;
	for (int i = 0; i < result.num; i++)
	{
		if (result.triElems[i].value > 1) result.triElems[i].value = 1;
	}
	return SUCCESS;
}

template<class T>
Status Relation_Symmetry(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& result)
{
	result.rows = tsma.rows; result.cols = tsma.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
	result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
	TriSparseMatrix<T> temp;tsma.FastTranspose(temp);
	Matrix_Add(tsma, temp, result);
	for (int i = 0; i < result.num; i++)
	{
		if (result.triElems[i].value > 1) result.triElems[i].value = 1;
	}
	return SUCCESS;
}

template<class T>
Status Warshall(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& result)
{
	result.rows = tsma.rows; result.cols = tsma.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
	result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
	T* p1 = tsma.Unseal();
	for (int i = 0; i < tsma.rows; i++)
	{
		T* A = new int[tsma.rows * tsma.cols];T* R = new int[tsma.rows]; T* C = new int[tsma.cols];
		for (int j = 0; j < tsma.cols; j++)
		{
			R[j] = *(p1 + i * tsma.cols + j);
		}
		for (int k = 0; k < tsma.rows; k++)
		{
			C[k] = *(p1 + k * tsma.cols + i);
		}
		Basic_Multiplication_numarray(C, R, tsma.rows, 1, 1, tsma.cols, A);
		Matrix_Add(A, p1, tsma.rows, tsma.cols, tsma.rows, tsma.cols, p1);
		delete[] R; delete[] C; delete[] A;
	}
	for (int i = 0; i < tsma.rows; i++)
	{
		for (int j = 0; j < tsma.cols; j++)
		{
			if (*(p1 + i * tsma.cols + j) != 0) result.SetElem(i, j, *(p1 + i * tsma.cols + j));
		}
	}
	for (int i = 0; i < result.num; i++)
	{
		if (result.triElems[i].value > 1) result.triElems[i].value = 1;
	}
	delete[] p1;
	return SUCCESS;
}

template <class T>
Status Matrix_Inversion(TriSparseMatrix<T>& tsma, TriSparseMatrix<T>& result)
{
	T* p1 = tsma.Unseal(); double determinant_value=numarray_determinant_value(p1,tsma.rows); 
	if (determinant_value == 0)
	{
		return FAIL;
	}
	double Algebraic_cofactor;
	T* inversion = new T[tsma.rows * tsma.cols]; T* sondet = new T[(tsma.rows - 1) * (tsma.cols - 1)];
	result.rows = tsma.rows; result.cols = tsma.cols; delete[] result.triElems; result.maxSize = DEFAULT_SIZE;
	result.triElems = new Triple<T>[result.maxSize]; result.num = 0;
	for (int i = 0; i < tsma.rows; i++)
	{
		for (int j = 0; j < tsma.cols; j++)
		{
			int ct = 0;
			for (int k = 0; k < tsma.rows; k++)
			{
				for (int l = 0; l < tsma.cols; l++)
				{
					if (k == i || l == j){continue;}
					else
					{
						*(sondet + ct) = *(p1 + k * tsma.cols + l); ct++;
					}
				}
			}
			Algebraic_cofactor = numarray_determinant_value(sondet, tsma.rows - 1);
			*(inversion + j * tsma.cols + i) = Algebraic_cofactor * pow(-1.0, (i + j + 2));
		}
	}
	for (int i = 0; i < tsma.rows; i++)
	{
		for (int j = 0; j < tsma.cols; j++)
		{

			*(inversion + i * tsma.cols + j)= *(inversion + i * tsma.cols + j)*(1.00/ determinant_value);
			if (*(inversion + i * tsma.cols + j) != 0)
			{
				result.SetElem(i,j, *(inversion + i * tsma.cols + j));
			}
		}
	}
	delete[] p1; delete[] sondet; delete[] inversion; return SUCCESS;
}

#endif
