// pca.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
using namespace std;


int main()
{
	const int feature_n = 2;
	const int cls_n = 2;
	int cls_sample_n[cls_n] = {4,3};
	const int total_n = 7;
	
	double data[feature_n][total_n] = {
	2.95, 2.53, 3.57, 3.16, 2.58, 2.16, 3.27,
	6.63, 7.79, 5.65, 5.47, 4.46, 6.22, 3.52
	};
	
	double** mydata = (double**) calloc(feature_n, sizeof(double*));
	for(int i = 0; i < feature_n; i++)
	{
		mydata[i] = (double*) calloc(total_n, sizeof(double));
	}
	for(int i = 0; i < feature_n; i++)
		for(int j = 0; j < total_n; j++)
		{
			mydata[i][j] = data[i][j];
		}
	
	double Eigen[feature_n*2]={0};
	double* EigenVector = (double*) calloc(feature_n*feature_n, sizeof(double));
	MutilClass_LDA(mydata, feature_n, cls_n, cls_sample_n, Eigen, EigenVector);

	system("pause");
}

//pca test
//int _tmain(int argc, _TCHAR* argv[])
//{
//	double data[20] = { 2, 3, 8, 2, 3,
//		7, 9, 29, 3, 5,
//		3, 8, 22, 12, 1,
//		3, 12, 12, 33, 2};
//	int row = 4;//特征维数
//	int col = 5;//样本个数
//	double signals[20], PC[16], V[4];
//	int x, y;
//	
//
//	pca1( data, row, col, signals, PC, V );
//	//pca2( data, row, col, signals, PC, V );
//
//	printf( "Project to Principal Component: \n" );
//	for ( y = 0; y < row; y++ )
//	{
//		for ( x = 0; x < col; x++ )
//		{
//			printf( "%7.4f ", signals[y*col+x] );
//		}
//		printf( "\n" );
//	}
//	printf( "Eigen Values (in deceasing order): \n" );
//	for ( y = 0; y < row; y++ )
//		printf( "%7.4f ", V[y] );
//	printf( "\n" );
//	printf( "Eigen Vectors: \n" );
//	for ( y = 0; y < row; y++ )
//	{
//		for ( x = 0; x < row; x++ )
//		{
//			printf( "%7.4f ", PC[y*row+x] );
//		}
//		printf( "\n" );
//	}
//
//
//	system("pause");
//	return 0;
//}
//
