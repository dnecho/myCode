// dianzhen.cpp : 定义控制台应用程序的入口点。

#include "opencv2/opencv.hpp"

#pragma comment(lib,"cvLib/opencv_highgui231.lib")
#pragma comment(lib,"cvLib/opencv_imgproc231.lib")
#include <stdio.h>
#include <tchar.h>
#include "dianzhen.h"

int _tmain(int argc, _TCHAR* argv[])
{	
	const char filename[255] = ".\\samples\\1.png";
	IplImage* pSrc = cvLoadImage(filename, 0);

	
	Axis center={0};
	int pixelsize=0;
	
	unsigned char* pSrcData = (unsigned char*)pSrc->imageData;
	
	unsigned char* pImgData = (unsigned char*)malloc(sizeof(unsigned char)*(pSrc->width*pSrc->height));
	for(int i = 0; i < pSrc->height; i++)
		for(int j = 0; j < pSrc->width; j++)
			pImgData[j+i*pSrc->width] = pSrcData[j+i*pSrc->widthStep];			
			
	process(pImgData, pSrc->width, pSrc->height, &center, &pixelsize);

	printf("final\t\tx:%d，y:%d，PixelSize:%d\n", center.x, center.y, pixelsize);
	
	system("pause");
	


	return 0;
}
