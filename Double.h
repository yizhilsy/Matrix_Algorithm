#include <iostream>
#include "Assistance.h"
using namespace std;
// ��Ԫ����
template<class ElemType>
class Double
{
public:
	// ���ݳ�Ա:
	int col;						// ����Ԫ�ص����±������±�
	ElemType value;						// ����Ԫ�ص�ֵ

	// ���캯��:
	Double();							// �޲����Ĺ��캯��
	Double(int c, ElemType v);	// ��֪����������Ԫ��
};

// ��Ԫ�����ʵ�ֲ���
template<class ElemType>
Double<ElemType>::Double()
// ����������������Ԫ��
{
}

template<class ElemType>
Double<ElemType>::Double(int c, ElemType v)
// �������������֪������������Ԫ��
{
	col = c;							// �к�
	value = v;							// ����Ԫ��ֵ
}
