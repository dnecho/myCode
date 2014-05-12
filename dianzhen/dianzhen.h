#pragma once

struct BaseCoordinate
{
	int pt1;
	int pt2;
	bool direction;
};

int PointInRegion(CPoint* center, int N, CPoint pt, double dis);

int FindBaseCoordinate(CPoint* center, int N, double distance, BaseCoordinate* baseCoordinate);

bool GetXyCoordinate(IplImage* pOutImg, int endPt1, int endPt2, CPoint* center, int N, double distance, bool direction);

void Detect(IplImage* pSrc, IplImage* pSrc3C);

void drawArrow(IplImage* img, CvPoint pStart, CvPoint pEnd, int len, int alpha,  CvScalar& color, int thickness=1, int lineType=4);