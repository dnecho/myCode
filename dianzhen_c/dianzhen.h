#ifndef DIANZHEN_H
#define DIANZHEN_H

#include "region.h"
#define REALDIS 0.33//��������Բ������ľ���Ϊ0.33����
typedef unsigned long DWORD;
#define des5_1Num 5//������5��1
#define ab_num 3
#define xynum 10//10����ʾ�����bitλ


#define MAX(a, b) (a)>(b)?(a):(b)
#define MIN(a, b) (a)<(b)?(a):(b)

typedef struct BaseCoordinate
{
	int pt1;
	int pt2;
	bool direction;
}BaseCoordinate;


typedef struct Axis {
	DWORD x;//x�᷽����ԭ��ľ��루����X������ֵ*1.98mm ����X�᷽��ƫ����������λ0.01mm��
	DWORD y;//y�᷽����ԭ��ľ��루����Y������ֵ*1.98mm ����Y�᷽��ƫ����������λ0.01mm��
} Axis;

void gaussianFilter(unsigned char* corrupted, unsigned char* smooth, int width, int height);

int PointInRegion(CPoint* center, int N, CPoint pt, double dis);

int FindBaseCoordinate(CPoint* center, int N, double distance, BaseCoordinate* baseCoordinate);

bool GetXyCoordinate(int endPt1, int endPt2, CPoint* PointCenter, int N, double distance, bool direction, CPoint PicCenter, double* tempX, double* tempY, double* PSize);

void adaptiveThreshold_C(unsigned char* input, int IMAGE_WIDTH, int IMAGE_HEIGHT, int IMAGE_WIDESTEP, unsigned char* bin);

void BilinearRowFilter(unsigned char* src, long* dst, int len, long* leftIdx, long* rightIdx, long* weight, int shift ,int nch);

void BiLinearInsert(unsigned char * pSrc, int sWidth, int sHeight, unsigned char * pDes, int dWidth, int dHeight ,int nch);

int process(unsigned char * yuy2buf, int width, int height, Axis *center, int *pixelsize);

#endif