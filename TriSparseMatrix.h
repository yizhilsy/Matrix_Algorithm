#ifndef __TRI_SPARSE_MATRIX_H__
#define __TRI_SPARSE_MATRIX_H__

#include "Assistance.h"				// 辅助软件包
#include "Triple.h"				// 三元组类
#include "Matrix_Function.h"//引入矩阵函数文件
template <typename ElemType>	class DouSparseMatrix;

// 稀疏矩阵三元组顺序表类
template<class ElemType>
class TriSparseMatrix
{
protected:
	// 稀疏矩阵三元组顺序表的数据成员:
	Triple<ElemType>* triElems;		// 存储稀疏矩阵的三元组表
	int maxSize;					// 非零元素最大个数
	int rows, cols, num;			// 稀疏矩阵的行数,列数及非零元个数

public:
	// 稀疏矩阵三元组顺序表的函数成员： 
	TriSparseMatrix(int rs = DEFAULT_SIZE, int cs = DEFAULT_SIZE, int size = DEFAULT_SIZE);// 构造一个rs行cs列非零元素最大个数为size的空稀疏矩阵
	TriSparseMatrix(ElemType* numarray,int rs,int cs , int size = DEFAULT_SIZE);	//从二维数组构建三元组表
	~TriSparseMatrix();				// 析构函数
	int GetRows() const;			// 返回稀疏矩阵行数
	int GetCols() const;			// 返回稀疏矩阵列数
	int GetNum() const;				// 返回稀疏矩阵非零元个数
	Status SetElem(int r, int c, const ElemType& v);	// 设置指定位置的元素值
	Status GetElem(int r, int c, ElemType& v);			// 求指定位置的元素值
	TriSparseMatrix(const TriSparseMatrix<ElemType>& copy);	// 复制构造函数
	TriSparseMatrix<ElemType>& operator =(const TriSparseMatrix<ElemType>& copy);// 赋值运算符重载
	void SimpleTranspose(TriSparseMatrix<ElemType>& b);// 稀疏矩阵的简单转置算法
	void FastTranspose(TriSparseMatrix<ElemType>& b);// 稀疏矩阵的快速转置算法
	ElemType* Unseal();//启封函数，将压缩的三元组在堆区解压后拷贝一份并返回地址
	void Print() const;//打印此三元组稀疏矩阵
	friend class DouSparseMatrix<ElemType>;//二元组类是三元组类的一个友类
	//矩阵/离散友元函数声明
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

// 稀疏矩阵三元组顺序表类的实现部分
template <class ElemType>
TriSparseMatrix<ElemType>::TriSparseMatrix(int r, int c, int size)
// 操作结果： 构造一个r行c列非零元素最大个数为size的空稀疏矩阵
{
	if (r < 1 || c < 1)
		throw Error("行数或列数无效!");	// 抛出异常
	maxSize = size;						// 非零元素最大个数
	rows = r;							// 行数
	cols = c;							// 列数
	num = 0;							// 非零元素个数
	triElems = new Triple<ElemType>[maxSize];	// 分配存储空间
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
// 操作结果：稀疏矩阵所占用空间
{
	if (triElems != NULL) delete[]triElems; // 释放存储空间
}

template <class ElemType>
int TriSparseMatrix<ElemType>::GetRows() const
// 操作结果：返回稀疏矩阵行数
{
	return rows;					// 返回行数
}

template <class ElemType>
int TriSparseMatrix<ElemType>::GetCols() const
// 操作结果：返回稀疏矩阵列数
{
	return cols;					// 返回列数
}

template <class ElemType>
int TriSparseMatrix<ElemType>::GetNum() const
// 操作结果：返回稀疏矩阵非零元个数
{
	return num;						// 返回非零元个数
}

template <class ElemType>
Status TriSparseMatrix<ElemType>::SetElem(int r, int c, const ElemType& v)
// 操作结果：如果下标范围错,则返回RANGE_ERROR,如果溢出,则返回OVER_FLOW,否则返
//	回SUCCESS
{
	if (r >= rows || c >= cols || r < 0 || c < 0)
		return RANGE_ERROR;					// 下标范围错

	int i, j;								// 工作变量
	for (j = num - 1; j >= 0 &&
		(r < triElems[j].row || r == triElems[j].row && c < triElems[j].col); j--);// 查找三元组位置

	if (j >= 0 && triElems[j].row == r && triElems[j].col == c) {	// 找到三元组
		if (v == 0) {	// 删除三元组
			for (i = j + 1; i < num; i++)
				triElems[i - 1] = triElems[i];	// 前移从j+1开始的三元组
			num--;								// 删除三元组后,非零元个数自减1
		}
		else	// 修改元素值
			triElems[j].value = v;
		return SUCCESS;						// 成功
	}
	else if (v != 0) {
		if (num < maxSize) {	// 将三元组(r, c, v)插入到三元组表中
			for (i = num - 1; i > j; i--)	// 后移元素	
				triElems[i + 1] = triElems[i];
			// j + 1为空出的插入位置
			triElems[j + 1].row = r;		// 行
			triElems[j + 1].col = c;		// 列
			triElems[j + 1].value = v;		// 非零元素值
			num++;							// 插入三元组后,非零元个数自加1
			return SUCCESS;					// 成功
		}
		else	// 溢出
			return OVER_FLOW;				// 溢出时返回OVER_FLOW
	}
	return SUCCESS;							// 成功
}

template <class ElemType>
Status TriSparseMatrix<ElemType>::GetElem(int r, int c, ElemType& v)
// 操作结果：如果下标范围错,则返回RANGE_ERROR,否则返回SUCCESS,并用v返回指定位置元素值
{
	if (r >= rows || c >= cols || r < 0 || c < 0)
		return RANGE_ERROR;			// 下标范围错

	int j;							// 工作变量


	for (j = num - 1; j >= 0 &&
		(r < triElems[j].row || r == triElems[j].row && c < triElems[j].col); j--);// 查找指定位置的三元组

	if (j >= 0 && triElems[j].row == r && triElems[j].col == c)	// 找到三元组
		v = triElems[j].value;		// 用v返回指定位置元素值
	else	// 未找到三元组
		v = 0;						// 未找到三元组,表示0元素值
	return SUCCESS;					// 成功
}

template <class ElemType>
TriSparseMatrix<ElemType>::TriSparseMatrix(const TriSparseMatrix<ElemType>& copy)
// 操作结果：由稀疏矩阵copy构造新稀疏矩阵——复制构造函数
{
	maxSize = copy.maxSize;							// 最大非零元素个数
	triElems = new Triple<ElemType>[maxSize];		// 分配存储空间
	rows = copy.rows;								// 复制行数
	cols = copy.cols;								// 复制列数
	num = copy.num;									// 复制非零元素个数
	triElems = new Triple<ElemType>[maxSize];		// 为三元组分配存储空间
	for (int i = 0; i < num; i++)	// 复制三元组
		triElems[i] = copy.triElems[i];
}

template <class ElemType>
TriSparseMatrix<ElemType>& TriSparseMatrix<ElemType>::operator =(const TriSparseMatrix<ElemType>& copy)
// 操作结果：将稀疏矩阵copy赋值给当前稀疏矩阵——赋值运算符重载
{
	if (&copy != this) {
		maxSize = copy.maxSize;						// 最大非零元素个数
		if (triElems != NULL) delete[]triElems;	// 释放存储空间
		triElems = new Triple<ElemType>[maxSize];	// 分配存储空间
		rows = copy.rows;							// 复制行数
		cols = copy.cols;							// 复制列数
		num = copy.num;								// 复制非零元素个数

		for (int i = 0; i < num; i++)	// 复制三元组
			triElems[i] = copy.triElems[i];
	}
	return *this;
}

template<class ElemType>
void TriSparseMatrix<ElemType>::SimpleTranspose(TriSparseMatrix<ElemType>& b)
// 操作结果：稀疏矩阵的简单转置算法，结果放在三元组顺序表b中 
{
	b.rows = cols;
	b.cols = rows;
	b.num = num;
	b.maxSize = maxSize;
	delete[]b.triElems;
	b.triElems = new Triple<ElemType>[b.maxSize];

	if (b.num > 0) {
		int i = 0;						// 稀疏矩阵b的下一个三元组的存放位置
		for (int col = 0; col < cols; col++)
			for (int j = 0; j < num; j++)	// 查找a矩阵中第col列的三元组
				if (triElems[j].col == col) {
					b.triElems[i].row = triElems[j].col;
					b.triElems[i].col = triElems[j].row;
					b.triElems[i].value = triElems[j].value;
					i++;
				}
	}
}

template<class ElemType>
void TriSparseMatrix<ElemType>::FastTranspose(TriSparseMatrix<ElemType>& b)//精妙算法！
// 操作结果：稀疏矩阵的快速转置算法，结果放在三元组顺序表b中 
{
	b.rows = cols;
	b.cols = rows;
	b.num = num;
	b.maxSize = maxSize;
	delete[]b.triElems;
	b.triElems = new Triple<ElemType>[b.maxSize];

	int* cNum = new int[cols + 1];	// 存放原矩阵中每一列的非零元个数
	int* cPos = new int[cols + 1];	// 存放原矩阵中每一列的第一个非零元在b中的存储位置
	int col;
	int i;

	if (b.num > 0) {
		for (col = 0; col < cols; col++) cNum[col] = 0;	// 初始化cNum
		for (i = 0; i < num; i++)
			++cNum[triElems[i].col];		// 统计原矩阵中每一列的非零元个数
		cPos[0] = 0;						// 第一列的第一个非零元在b存储的起始位置
		for (col = 1; col < cols; col++)	// 循环求每一列的第一个非零元在b存储的起始位置
			cPos[col] = cPos[col - 1] + cNum[col - 1];

		for (i = 0; i < num; i++) {	// 循环遍历原矩阵的三元组
			int j = cPos[triElems[i].col];
			// 用于表示b当前列的下一个非零元三元组的存储位置
			b.triElems[j].row = triElems[i].col;
			b.triElems[j].col = triElems[i].row;
			b.triElems[j].value = triElems[i].value;
			++cPos[triElems[i].col];
			// b当前列的下一个非零元三元组的存储新位置
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
	return point;//返回在堆区开辟空间指针，在函数外使用后需要记得删除
}

template<class ElemType>
void TriSparseMatrix<ElemType>::Print() const
{
	cout << "三元组表现形式:" << endl;
	for (int i = 0; i < num; i++)
	{
		cout << "(" << triElems[i].row << "," << triElems[i].col << "," << triElems[i].value << ");";
	}
	cout << endl << "矩阵表现形式:" << endl; int flag = 0;
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

//矩阵/离散友元函数定义
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
		delete[] p1; delete[] p2; delete[] presult;//删除在堆区解压的元素 
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
		delete[] p1; delete[] p2; delete[] presult;//删除在堆区解压的元素 
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
		//计算row_offsets以及r_num辅助计算矩阵的乘积
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
		//利用row_offsets以及r_num计算矩阵的乘积
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
		delete[] p2; //删除在堆区解压的元素 
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
		delete[] p1; delete[] p2; delete[] presult;//删除在堆区解压的元素 
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
		delete[] p1; delete[] p2; delete[] presult;//删除在堆区解压的元素 
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
