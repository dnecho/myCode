#ifndef REGION_H
#define REGION_H

#define bool int
#define false 0
#define true 1

typedef struct Region
{
	int PointNum;
	int left;
	int right;
	int top;
	int bottom;
	int FusedN;
	int W;
	int H;

	int LeftTop;
	int RightTop;

	bool IsOK;
}Region;

typedef struct CPoint
{
	int x;
	int y;
}CPoint;

typedef struct RectCenter{
	double x;
	double y;
}RectCenter;

int KeyGrow(unsigned char * p, int w, int h,int * typeImg);

#endif