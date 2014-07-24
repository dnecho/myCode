// fapiao0721.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "region.h"
IplImage* filterSigleChannel(IplImage* pSrc, char chn, int RTh, int GTh, int BTh);
IplImage* detectEllipse(IplImage* pSrc);

int _tmain(int argc, _TCHAR* argv[])
{

	IplImage* pImg = cvLoadImage("samples//图片3.jpg");

	//IplImage* pGrey = cvCreateImage(cvGetSize(pImg), pImg->depth, 1);
	//cvCvtColor(pImg, pGrey, CV_RGB2GRAY);
	if(!pImg)
	{
		printf("open image error");
		system("pause");
		return 0;
	}

	IplImage* pWinImg = detectEllipse(pImg);
	cvSmooth(pWinImg, pWinImg);
	IplImage* pTemp = filterSigleChannel(pWinImg, 'r', 130, 180, 180);



	cvNamedWindow("src");
	cvShowImage("src", pImg);
	cvSaveImage("out.jpg",pWinImg);
	cvNamedWindow("out");
	cvShowImage("out", pTemp);
	cvWaitKey();

	//system("pause");
	return 0;
}

IplImage* filterSigleChannel(IplImage* pSrc, char chn, int RTh, int GTh, int BTh)
{
	IplImage *pSrcR = cvCreateImage(cvGetSize(pSrc), pSrc->depth, 1);
	IplImage *pSrcG = cvCreateImage(cvGetSize(pSrc), pSrc->depth, 1);
	IplImage *pSrcB = cvCreateImage(cvGetSize(pSrc), pSrc->depth, 1);
	cvSplit(pSrc, pSrcB, pSrcG, pSrcR, NULL);

	unsigned char* pData = (unsigned char*)pSrc->imageData;
	unsigned char rgb[3];
	int th[3];
	IplImage *pTemp = cvCreateImage(cvGetSize(pSrc), pSrc->depth, 1);
	unsigned char* pTempData = (unsigned char*)pTemp->imageData;

	for(int i = 0; i < pSrc->height; i++)
	{
		for(int j = 0; j < pSrc->width; j++)
		{
			int pos = i*pSrc->widthStep+j*3;
			switch(chn){
				case 'r':
					th[2]=RTh; th[1]=GTh; th[0]=BTh;
					rgb[2] = pData[pos+2];//R
					rgb[1] = pData[pos+1];//G
					rgb[0] = pData[pos];//B
					break;
				case 'g':
					th[2]=GTh; th[1]=BTh; th[0]=RTh;
					rgb[2] = pData[pos+1];//G
					rgb[1] = pData[pos];//B
					rgb[0] = pData[pos+2];//R
					break;
				case 'b':
					th[2]=BTh; th[1]=GTh; th[0]=RTh;
					rgb[2] = pData[pos];//B
					rgb[1] = pData[pos+1];//G
					rgb[0] = pData[pos+2];//R
					break;
				default:
					cvReleaseImage(&pSrcR);
					cvReleaseImage(&pSrcG);
					cvReleaseImage(&pSrcB);
					return pTemp;

			}

			//if( rgb[1]<gbT && rgb[0]<gbT && rgb[2]>rT)
			if( rgb[0]<th[0] && rgb[1]<th[1] && rgb[2]>th[2])
			{
				pTempData[i*(pTemp->widthStep)+j] = 0;
			}else
			{
				pTempData[i*(pTemp->widthStep)+j] = 255;
			}
		}
	}
	cvReleaseImage(&pSrcR);
	cvReleaseImage(&pSrcG);
	cvReleaseImage(&pSrcB);
	return pTemp;
}

IplImage* detectEllipse(IplImage* pSrc)
{
	int width = 500;
	float ratio = 1.0*pSrc->width/width;
	int height = (int)(pSrc->height/ratio);

	IplImage* pSrc3C = cvCreateImage(cvSize(width,height), pSrc->depth, pSrc->nChannels);
	cvResize(pSrc, pSrc3C);
	cvSmooth(pSrc3C, pSrc3C, CV_GAUSSIAN);
	int widthStep = pSrc3C->widthStep;


	//1.过滤红色通道的图章信息，得到2值图像
	IplImage* pTemp = filterSigleChannel(pSrc3C, 'r', 150, 150, 150);
	unsigned char* pTempData = (unsigned char*)pTemp->imageData;

	//2.去除小区域，避免干扰
	int* typeImg=(int*)malloc(sizeof(int)*height*width);

	extern Region region[REGION_NUM];
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			typeImg[i*width+j]=-1;
		}
	}
	int regionN = KeyGrow(pTempData, width, height, typeImg);
	cout<<"检测出的连通区域个数:排除小区域后剩余区域个数\t\t"<<regionN<<":";
	pTempData = (unsigned char*)pTemp->imageData;
	int tempN=0;
	for(int i = 0; i < regionN; i++)
	{
		if(region[i].PointNum < 500)
		{
			tempN++;
			for(int ih=region[i].top; ih <=region[i].bottom; ih++)
			{
				for(int iw=region[i].left; iw<=region[i].right; iw++)
				{
					pTempData[iw+ih*pTemp->widthStep]=255;
				}
			}
		}
	}
	cout<<tempN<<endl;

	//3.检查椭圆
	CvMemStorage* storage = cvCreateMemStorage(0);  
	CvSeq* contours;  
	CvBox2D s;  
	int num = cvFindContours( pTemp, storage, &contours, sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cvPoint(0,0));  
	//s=cvFitEllipse2(contours);
	int ElliNum = 0;


	CvBox2D bottomBox={0};
	cout<<"\nangel\tx\ty\tw\th"<<endl;
	for(CvSeq* c=contours; c!=NULL; c=c->h_next)  
	{  
		if(c->total<150)
			continue;
		s=cvFitEllipse2(c);

		if(s.center.y>bottomBox.center.y)
		{
			bottomBox=s;
		}
		ElliNum++;

		cout<<s.angle<<"\t"<<s.center.x<<"\t"<<s.center.y<<"\t"<<s.size.width<<"\t"<<s.size.height<<endl;

	}
	cout<<"\n检测出的椭圆个数：排除小椭圆后剩余椭圆个数\t\t"<<num<<":"<<ElliNum<<endl;
	if(!ElliNum)
	{
		cout<<"can't detect ellipse"<<endl;
		return NULL;
	}
	//cvEllipseBox(pSrc3C, bottomBox, CV_RGB(0,0,255),2, 8 , 0 );//画椭圆框

	//4.根据检测到的椭圆映射回原有图像（旋转后），并提取图章返回Img
	bottomBox.angle = bottomBox.angle>180 ? bottomBox.angle-180 : bottomBox.angle;
	bottomBox.center.x *= ratio;
	bottomBox.center.y *= ratio;
	bottomBox.size.height *= ratio;
	bottomBox.size.width *= ratio;

	CvPoint2D32f pt[4];
	CvPoint2D32f RoPt[4];
	cvBoxPoints(bottomBox,pt);

	float angle = bottomBox.angle-90;
	float sina = sin(RADIAN(angle));
	float cosa = cos(RADIAN(angle));
	IplImage* pRoImg = RotateImage_CV(pSrc, angle);//顺时间方向转动为负角度
	for(int i=0; i<4; ++i)
	{
		RoPt[i].x = (pt[i].x-pSrc->width/2)*cosa + (pt[i].y-pSrc->height/2)*sina + pRoImg->width/2;
		RoPt[i].y = -(pt[i].x-pSrc->width/2)*sina + (pt[i].y-pSrc->height/2)*cosa + pRoImg->height/2;
	}

	int rect[4];//rect:x,y,w,h
	memset(rect, 1,sizeof(int)*4);
	CvPoint rectPt[4]={0};
	for(int i = 0; i < 4; i++)
	{
		rectPt[i] = cvPointFrom32f(RoPt[i]);
		rect[0] = (rect[0]>rectPt[i].x) ? rectPt[i].x : rect[0];//x
		rect[1] = (rect[1]>rectPt[i].y) ? rectPt[i].y : rect[1];//y
	}
	rect[2]=0;
	rect[3]=0;
	for(int i = 0; i < 3; i++)
	{
		int rectw = abs(rectPt[i].x - rectPt[i+1].x);
		int recth = abs(rectPt[i].y - rectPt[i+1].y);
		rect[2] = (rectw>rect[2]) ? rectw : rect[2];//w
		rect[3] = (recth>rect[3]) ? recth : rect[3];//h
	}

	CvRect winRect = cvRect(rect[0], rect[1], rect[2], rect[3]);
	cvSetImageROI(pRoImg, winRect);
	IplImage* pWinImg = cvCreateImage(cvSize(winRect.width, winRect.height), pRoImg->depth, pRoImg->nChannels);
	cvCopy(pRoImg, pWinImg);
	//IplImage* pWinImg = cvCloneImage(pRoImg);
	cvResetImageROI(pRoImg);

	//cvNamedWindow("pSrc3C");
	//cvShowImage("pSrc3C", pWinImg);

	cvReleaseImage(&pSrc3C);
	cvReleaseMemStorage(&storage);
	delete []typeImg;
	return pWinImg;
}
