#pragma once

//////////////////////////////////////////////////////////////////////////
//函数用途：多类别LDA，求得特征值和特征向量
//in:
//		data:样本数据，按类别顺序依次排列，特征为列方向，样本为行方向
//		feature_n:特征种类
//		cls_n:类别数
//		cls_sample_n:每个类别的样本个数
//out:
//		Eigen:求得的特征值,在函数外申请内存，大小为feature_n
//		EigenVector:求得的特征向量，在函数外申请内存，大小为feature_n*feature_n
//////////////////////////////////////////////////////////////////////////
void MutilClass_LDA(double** data,int feature_n, int cls_n, int* cls_sample_n, double* Eigen, double* EigenVector);

//////////////////////////////////////////////////////////////////////////
//函数用途：把方阵化为上三角Hessenberg矩阵
//A1:n阶方阵
//
//n:方阵A1的阶数
//
//Ret:返回的n阶方阵
//////////////////////////////////////////////////////////////////////////
void Matrix_Hessenberg(double *A1,int n,double *ret);

//////////////////////////////////////////////////////////////////////////
//函数用途：通过QR分解法求矩阵特征值
//K1：n阶方阵
//n：方阵K1的阶数
//LoopNumber：在误差无法保证能得到结果时运算的最大次数
//Error1：误差控制变量
//Ret：返回的一个n*2的矩阵。矩阵的每一行是求得的特征值，第一列代表特征值实数部分，第二列代表虚数部分
//函数成功返回True，失败返回False
//////////////////////////////////////////////////////////////////////////
bool Matrix_EigenValue(double *K1,int n,int LoopNumber,double Error1,double *Ret);