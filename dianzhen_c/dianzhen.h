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
	DWORD x;//x轴方向离原点的距离（等于X轴坐标值*1.98mm 加上X轴方向偏移量），单位0.01mm，
	DWORD y;//y轴方向离原点的距离（等于Y轴坐标值*1.98mm 加上Y轴方向偏移量），单位0.01mm，
} Axis;

void gaussianFilter(unsigned char* corrupted, unsigned char* smooth, int width, int height);

int PointInRegion(CPoint* center, int N, CPoint pt, double dis);

int FindBaseCoordinate(CPoint* center, int N, double distance, BaseCoordinate* baseCoordinate);

bool GetXyCoordinate(int endPt1, int endPt2, CPoint* PointCenter, int N, double distance, bool direction, CPoint PicCenter, double* tempX, double* tempY, double* PSize);

void adaptiveThreshold_C(BYTE* input, int IMAGE_WIDTH, int IMAGE_HEIGHT, int IMAGE_WIDESTEP, BYTE* bin);

int process(unsigned char * yuy2buf, int width, int height, Axis *center, int *pixelsize);

