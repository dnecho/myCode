#pragma once
typedef  unsigned char BYTE;

struct Region
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
};

struct CPoint
{
	int x;
	int y;
};

struct RectCenter{
	double x;
	double y;
};

int KeyGrow(unsigned char * p, int w, int h,int * typeImg);