#include <iostream>
#include "Assistance.h"
using namespace std;
// 三元组类
template<class ElemType>
class Double
{
public:
	// 数据成员:
	int col;						// 非零元素的行下标与列下标
	ElemType value;						// 非零元素的值

	// 构造函数:
	Double();							// 无参数的构造函数
	Double(int c, ElemType v);	// 已知数据域建立三元组
};

// 三元组类的实现部分
template<class ElemType>
Double<ElemType>::Double()
// 操作结果：构造空三元组
{
}

template<class ElemType>
Double<ElemType>::Double(int c, ElemType v)
// 操作结果：由已知数数据域构造三元组
{
	col = c;							// 列号
	value = v;							// 非零元素值
}
