#include "region.h"
#define REALDIS 0.33
typedef unsigned long DWORD;

struct BaseCoordinate
{
	int pt1;
	int pt2;
	bool direction;
};


typedef struct Axis {
	DWORD x;//x�᷽����ԭ��ľ��루����X������ֵ*1.98mm ����X�᷽��ƫ����������λ0.01mm��
	DWORD y;//y�᷽����ԭ��ľ��루����Y������ֵ*1.98mm ����Y�᷽��ƫ����������λ0.01mm��
} Axis;

void gaussianFilter(unsigned char* corrupted, unsigned char* smooth, int width, int height);

int PointInRegion(CPoint* center, int N, CPoint pt, double dis);

int FindBaseCoordinate(CPoint* center, int N, double distance, BaseCoordinate* baseCoordinate);

bool GetXyCoordinate(int endPt1, int endPt2, CPoint* PointCenter, int N, double distance, bool direction, CPoint PicCenter, double* tempX, double* tempY, double* PSize);

void adaptiveThreshold_C(BYTE* input, int IMAGE_WIDTH, int IMAGE_HEIGHT, int IMAGE_WIDESTEP, BYTE* bin);

int process(unsigned char * yuy2buf, int width, int height, Axis *center, int *pixelsize);

