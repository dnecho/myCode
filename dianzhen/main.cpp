// dianzhen.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"

int _tmain(int argc, _TCHAR* argv[])
{	
	const char filename[255] = ".\\samples\\Image011.bmp";
	IplImage* pSrc = cvLoadImage(filename, 0);
	IplImage* pSrc3C = cvLoadImage(filename, 1);
	Detect(pSrc, pSrc3C);

	return 0;
}
