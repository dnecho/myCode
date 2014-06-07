#include "stdafx.h"
#include "dianzhen.h"
#include "region.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <memory.h>

int g_ImgX = 0;
int g_ImgY = 0;
 bool IsPointInImg(CPoint pt)
{
	if(pt.x>=0 && pt.y>=0 && pt.x<g_ImgX && pt.y<g_ImgY)
		return true;
	return false;
}

 double Get2PointDis(CPoint a, CPoint b)
{
	return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
}

int PointInRegion(CPoint* center, int N, CPoint pt, double dis)
{
	double d = dis*0.5;
	if(!IsPointInImg(pt))//����õ㲻��ͼ��
		return -2;
	for(int i = 0; i < N; i++)
	{
		double temp = Get2PointDis(pt,center[i]);
		if(temp < d)
			return i;
	}
	return -1;
}


int FindBaseCoordinate(CPoint* center, int N, double distance, BaseCoordinate* baseCoordinate)
{
	int BaseCoordinateNum = 0;
	int endpt1_idx, endpt2_idx;
	double delta_4x,delta_4y;

	double sigma = 1.8;//��Ӧ�Ȳ������ò���Ӧ��С��2
	const int des5_1Num = 5;//��5��1
	int pt5Idx[des5_1Num]={0};
	for(int i = 0; i < N; i++)
	{
		for(int j = i+1; j < N; j++)
		{
			//if(i == 22)
			//	int mydebug = 1;
			//int endpt1_idx,endpt2_idx;//�ҵ�5������1�ĵ��е������˵�
			int curridx,nextInResion;
			double tempDis = Get2PointDis(center[i], center[j]);
			if(tempDis > sigma*distance)//�˵���Χû�������㣬����
				continue;			

			endpt1_idx = i;
			endpt2_idx = j;

			int delta_x = center[j].x - center[i].x;
			int delta_y = center[j].y - center[i].y;
			pt5Idx[0] = i;
			pt5Idx[1] = j;
			int ptNum = 0;

			//��������Ŀ���1
			CPoint beforePt;
			beforePt.x = center[i].x - delta_x;
			beforePt.y = center[i].y - delta_y;
			curridx = i;
			nextInResion = PointInRegion(center, N, beforePt,distance);
			while(nextInResion>=0 && ptNum < 3)//��һ����ĵ���Ŀ���1
			{
				endpt1_idx = nextInResion;
				pt5Idx[ptNum+2] = nextInResion;
				ptNum++;
				beforePt.x = 2*center[nextInResion].x - center[curridx].x;
				beforePt.y = 2*center[nextInResion].y - center[curridx].y;

				curridx = nextInResion;
				nextInResion = PointInRegion(center, N, beforePt,distance);
			}

			//�����Ҳ��Ŀ���1
			CPoint nextPt;
			nextPt.x = center[j].x+delta_x;
			nextPt.y = center[j].y+delta_y;
			curridx = j;
			nextInResion = PointInRegion(center, N, nextPt,distance);
			while(nextInResion>=0 && ptNum < des5_1Num-2)//��һ����ĵ���Ŀ���1
			{
				endpt2_idx = nextInResion;
				pt5Idx[ptNum+2] = nextInResion;
				ptNum++;
				nextPt.x = 2*center[nextInResion].x - center[curridx].x;
				nextPt.y = 2*center[nextInResion].y - center[curridx].y;

				curridx = nextInResion;
				nextInResion = PointInRegion(center, N, nextPt,distance);
			}

			if(ptNum!=des5_1Num-2)//����������5��������1�ſ���������
				continue;

			//�ж϶˵��Ƿ��Ѿ����ҳ���
			bool bFindBaseCoordinate = false;
			for(int BaseCoordinate_i = 0; BaseCoordinate_i < BaseCoordinateNum; BaseCoordinate_i++)
			{
				if(baseCoordinate[BaseCoordinate_i].pt1 == endpt1_idx || baseCoordinate[BaseCoordinate_i].pt1 == endpt2_idx)
				{
					bFindBaseCoordinate = true;
					break;
				}
			}
			if(bFindBaseCoordinate)
				continue;


			//�ж������˵������Ƿ�Ϊ0������������5���㲻������
			delta_4x = (center[endpt2_idx].x - center[endpt1_idx].x)/4.0;
			delta_4y = (center[endpt2_idx].y - center[endpt1_idx].y)/4.0;
			CPoint pt1;
			pt1.x = (int)(center[endpt1_idx].x-delta_4x);
			pt1.y = (int)(center[endpt1_idx].y-delta_4y);
			if(PointInRegion(center, N, pt1, distance) != -1)//�ж��Ƿ�Ϊ0������������5���㲻������,����
				continue;
			CPoint pt2;
			pt2.x = (int)(center[endpt2_idx].x+delta_4x);
			pt2.y = (int)(center[endpt2_idx].y+delta_4y);
			int tep = PointInRegion(center, N, pt2, distance);
			if(PointInRegion(center, N, pt2, distance) != -1)//�ж��Ƿ�Ϊ0������������5���㲻������,����
				continue;

			//��5��������pt1,pt2��0�Ļ����ϣ��ж�pt1,pt2�����Ƿ�Ϊ10/11
			const int ab_num = 2;
			CPoint a[ab_num],b[ab_num];//(above,below)
			for(int ab_i = 0; ab_i < ab_num; ab_i++)
			{
				a[ab_i].x = (int)(pt1.x + delta_4y);
				a[ab_i].y = (int)(pt1.y - delta_4x);

				b[ab_i].x = (int)(pt1.x - delta_4y);
				b[ab_i].y = (int)(pt1.y + delta_4x);
			}

			if(PointInRegion(center, N, a[1], distance)<0 && PointInRegion(center, N, b[1], distance)<0)//��a1,b1���붼Ϊ1����������
				continue;
			int a0_Region = PointInRegion(center, N, a[0], distance);
			int b0_Region = PointInRegion(center, N, b[0], distance);

			if(a0_Region==-1 && b0_Region>=0 )
			{
				baseCoordinate[BaseCoordinateNum].pt1 = endpt1_idx;
				baseCoordinate[BaseCoordinateNum].pt2 = endpt2_idx;
				baseCoordinate[BaseCoordinateNum].direction = true;
				BaseCoordinateNum++;
			}
			else if(a0_Region>0 && b0_Region==-1 )
			{
				baseCoordinate[BaseCoordinateNum].pt1 = endpt1_idx;
				baseCoordinate[BaseCoordinateNum].pt2 = endpt2_idx;
				baseCoordinate[BaseCoordinateNum].direction = false;
				BaseCoordinateNum++;
			}				
		}
	}

	return BaseCoordinateNum;
}

//���������ڻ�
double DotProduct(CPoint a, CPoint b)
{
	return a.x*b.x + a.y*b.y;
}


bool GetXyCoordinate(int endPt1, int endPt2, CPoint* PointCenter, int N, double distance, bool direction, CPoint PicCenter, double* tempX, double* tempY, double* PSize)
{
	CPoint BasePt = PointCenter[endPt1];
	double delta_x = (PointCenter[endPt2].x - PointCenter[endPt1].x)/4.0;
	double delta_y = (PointCenter[endPt2].y - PointCenter[endPt1].y)/4.0;

	const int xynum = 10;
	CPoint xPt[xynum],yPt[xynum];
	int x_ratio,y_ratio;//����ϵ����ͬ
	for(int i = 0; i < xynum; i++)
	{	
		x_ratio = i>4?1:2;
		xPt[i].x = (int)(BasePt.x + x_ratio*delta_y + (i%5)*delta_x);
		xPt[i].y = (int)(BasePt.y - x_ratio*delta_x + (i%5)*delta_y);

		y_ratio = i>4?2:1;
		yPt[i].x = (int)(BasePt.x - y_ratio*delta_y + (i%5)*delta_x);
		yPt[i].y = (int)(BasePt.y + y_ratio*delta_x + (i%5)*delta_y);

	}

	//����ĸ��˵㲻��ͼ���У�������
	if(!IsPointInImg(xPt[0]) || !IsPointInImg(xPt[4]) || !IsPointInImg(yPt[5]) || !IsPointInImg(yPt[9]))
		return false;

	int xOutAarry[xynum],yOutAarry[xynum];
	for(int i = 0; i < xynum; i++)
	{
		xOutAarry[i] = PointInRegion(PointCenter, N, xPt[i], distance)>=0?1:0;
		yOutAarry[i] = PointInRegion(PointCenter, N, yPt[i], distance)>=0?1:0;
	}

	if(!direction)//������򲻶ԣ������������
	{
		int xytmp;
		for(int i = 0; i < xynum; i++)
		{
			xytmp = xOutAarry[i];
			xOutAarry[i] = yOutAarry[9-i];
			yOutAarry[9-i] = xytmp; 
		}
	}

	int xOut = 0;
	int yOut = 0;

	for(int i = 0; i < xynum; i++)
	{
		xOut = (xOut<<1) | xOutAarry[i];
		yOut = (yOut<<1) | yOutAarry[i];
	}
	printf("����x:%d,y:%d\t\t",xOut, yOut );

	//����ԭʼ����ϵx��y�����᷽��
	CPoint StartPt,EndPt;
	
	StartPt.x = direction?xPt[7].x:yPt[2].x;
	StartPt.y = direction?xPt[7].y:yPt[2].y;
	EndPt.x = direction?xPt[9].x:yPt[0].x;
	EndPt.y = direction?xPt[9].y:yPt[0].y;

	CPoint V_X;//={(), (EndPt.y-StartPt.y)};//ԭʼx��������
	V_X.x=EndPt.x-StartPt.x;
	V_X.y=EndPt.y-StartPt.y;
	
	StartPt.x = direction?xPt[7].x:yPt[2].x;
	StartPt.y = direction?xPt[7].y:yPt[2].y;
	EndPt.x = direction?yPt[2].x:xPt[7].x;
	EndPt.y = direction?yPt[2].y:xPt[7].y;

		
	CPoint V_Y;//={(EndPt.x-StartPt.x), (EndPt.y-StartPt.y)};//ԭʼy��������
	V_Y.x=EndPt.x-StartPt.x;
	V_Y.y=EndPt.y-StartPt.y;
	
	double xyDotP = DotProduct(V_X, V_Y);//��xy�����������ڻ�����Ϊx,y���괹ֱ������ֵӦ��Ϊ0��������ܷ���˵����������ȷ��
	
	//��xy����������ģ
	double V_Xmode = sqrt(1.0*(V_X.x*V_X.x + V_X.y*V_X.y));
	double V_Ymode = sqrt(1.0*(V_Y.x*V_Y.x + V_Y.y*V_Y.y));
		

	//��ⵥԪ�����ĵ�
	CPoint UnitCenter;
	UnitCenter.x = (xPt[7].x + yPt[2].x)/2;
	UnitCenter.y = (xPt[7].y + yPt[2].y)/2;
	

	CPoint V_U2P;//��ⵥԪ�����ĵ���ͼƬ������ɵ��������ɼ�ⵥԪ������ָ��ͼƬ����
	V_U2P.x = PicCenter.x - UnitCenter.x;
	V_U2P.y = PicCenter.y - UnitCenter.y;
	

	
	//��������V_U2P��x��y�����ͶӰ
	double x_shadow = DotProduct(V_X, V_U2P)/V_Xmode;
	double y_shadow = DotProduct(V_Y, V_U2P)/V_Ymode;
	
	double PixelsDis_2Point = sqrt(delta_x*delta_x + delta_y*delta_y);
	(*PSize) = REALDIS/PixelsDis_2Point ;//���ش����ʵ�ʾ���
	
	(*tempX) = (*PSize) * x_shadow + xOut*1.98;//x�����ƫ����
	(*tempY) = (*PSize) * y_shadow + yOut*1.98;//y�����ƫ����
	
	return true;
}

int small_gaussian[4][4] =
{
	{1024}, 
	{512, 256}, 
	{384 , 256, 64},
	{288, 224, 112, 32}
};


//��ͨ��ͼ���˹ƽ��
void MyGaussian(unsigned char *pSrcData, unsigned char *pDesData, int Width, int Height)
{
	int i,j,k,m;
	//int Width=pSrc->width;
	//int Height=pSrc->height;
	int nCh=1;
	int widestep=Width;
	BYTE * Data=pSrcData;	
	BYTE * sData;
	int SmoothRadio=3;
	int pos,rpos;
	int Temp;
	int *  pGauTepData=NULL;

	//����SMOOTHRADIO�ж����ĸ�����
	if (SmoothRadio==1)
	{
		pGauTepData=small_gaussian[1];
	}
	else if(SmoothRadio==2)
	{
		pGauTepData=small_gaussian[2];
	}
	else if(SmoothRadio==3)
	{
		pGauTepData=small_gaussian[3];
	}
	else
	{

	}

	BYTE * rfilter=(BYTE*)malloc(sizeof(BYTE)*Width*Height);
	sData=rfilter;
	//���˲�
	pos=0;
	rpos=0;
	for (i=0;i<Height;i++)
	{
		for (j=3;j<Width-3;j++)
		{
			rpos=pos+j*nCh;
			for (m=0;m<nCh;m++)
			{
				Temp=Data[rpos+m]*pGauTepData[0];
				for (k=1;k<=SmoothRadio;k++)
				{
					Temp+=(Data[rpos+m-k*nCh]+Data[rpos+m+k*nCh])*pGauTepData[k];
				}

				sData[rpos+m]=Temp>>10;
			}
		}
		pos+=widestep;
	}

	sData=pDesData;
	Data=rfilter;
	pos=3*widestep;
	rpos=0;
	for (i=3;i<Height-3;i++)
	{
		for (j=3;j<Width-3;j++)
		{
			rpos=pos+j*nCh;
			for (m=0;m<nCh;m++)
			{
				Temp=Data[rpos+m]*pGauTepData[0];
				for (k=1;k<=SmoothRadio;k++)
				{
					Temp+=(Data[rpos+m-k*widestep]+Data[rpos+m+k*widestep])*pGauTepData[k];
				}

				sData[rpos+m]=Temp>>10;
			}
		}
		pos+=widestep;
	}

	free(rfilter);
}

void cvText(IplImage* img, const char* text, int x, int y)
{
	CvFont font;

	double hscale = 0.3;
	double vscale = 0.3;
	int linewidth = 1;
	cvInitFont(&font,CV_FONT_HERSHEY_SIMPLEX ,hscale,vscale,0,linewidth);

	CvScalar textColor =cvScalar(0,255,255);
	CvPoint textPos =cvPoint(x, y);

	cvPutText(img, text, textPos, &font,textColor);
}

int process(unsigned char * yuy2buf, int width, int height, Axis *center, int *pixelsize, CString& Outcome)
{

	g_ImgX = width;
	g_ImgY = height;

	//1.��˹�˲�
	unsigned char* pDataGaussian = (unsigned char*)malloc(sizeof(unsigned char)*(width*height));

	memcpy(pDataGaussian, yuy2buf, sizeof(unsigned char)*(width*height));

	//gaussianFilter(yuy2buf, pDataGaussian, width, height);
	MyGaussian(yuy2buf, pDataGaussian, width, height);

	//2.��ֵ��
	unsigned char* pData = (unsigned char*)malloc(sizeof(BYTE)*width*height);
	adaptiveThreshold_C(pDataGaussian, width, height, width, pData);
	

	//3.Ѱ����ͨ����
	int * typeImg=(int*)malloc(sizeof(int)*height*width);
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			typeImg[i*width+j]=-1;
		}
	}
	int N=KeyGrow(pData,width,height,typeImg);

	extern Region region[1000];
	char txt[255]="";

	//4.ͨ������С���������ȷ��դ�����distance
	//IplImage* pShow = cvLoadImage(".//samples//Image03.bmp",1);//FOR_TEST
	//char NumChar[255];//FOR_TEST
	CPoint* PointCenter = (CPoint*)malloc(sizeof(CPoint)*N);//����ͨ�����С��Ӿ��������
	for(int i = 0; i < N; i++)
	{
		PointCenter[i].x = (region[i].left+region[i].right)>>1;
		PointCenter[i].y = (region[i].top+region[i].bottom)>>1;
		//if(i==8||i==19)
		//{
		//sprintf_s(NumChar,255, "%d",i);
		//cvText(pShow,NumChar, PointCenter[i].x, PointCenter[i].y);
		//cvRectangle(pShow, cvPoint(region[i].left, region[i].top), cvPoint(region[i].right, region[i].bottom), cvScalar(255,0,255));//FOR_TEST
		//}
	}
	
	//////////////////////////////////////////////////////////////////////////
	//cvNamedWindow("test",0);
	//cvShowImage("test", pShow);
	//cvWaitKey();
	//////////////////////////////////////////////////////////////////////////
	
	
	double distance = 0xFFFFFFF;
	double boxRatio = 0.2;//����̾���ʱȥ����ͼƬ���ܵĵ㣬��ֹ�ҵ������diastance
	int leftBox = boxRatio * width;
	int rightBox = width - leftBox;
	int topBox = boxRatio * height;
	int bottomBox = height - topBox;
	for(int i = 0; i < N; i++)
	{
		if(PointCenter[i].x < leftBox || PointCenter[i].x > rightBox || PointCenter[i].y < topBox || PointCenter[i].y > bottomBox)
			continue;
		for(int j = i+1; j < N; j++)
		{
			if(PointCenter[j].x < leftBox || PointCenter[j].x > rightBox || PointCenter[j].y < topBox || PointCenter[j].y > bottomBox)
				continue;
			double tempDis = Get2PointDis(PointCenter[i], PointCenter[j]);
			if(tempDis < distance)
			{
				distance = tempDis;
			}
		}
	}

	CPoint PicCenter;//={width/2, height/2};
	PicCenter.x=width/2;
	PicCenter.y=height/2;

	//5Ѱ�һ�����
	BaseCoordinate* baseCoordinate = (BaseCoordinate*)malloc(sizeof(BaseCoordinate)*N/8);//ÿ��դ��������8��1
	int BaseCoordinateNum = FindBaseCoordinate(PointCenter, N, distance, baseCoordinate);
	int SuccessUnitNum = 0;

	center->x = 0;
	center->y = 0;
	double CenterX = 0;
	double CenterY = 0;
	double Psize = 0;
	double tempX,tempY,tempPsize;
	for(int i = 0; i < BaseCoordinateNum; i++)
	{
		bool success = GetXyCoordinate(baseCoordinate[i].pt1, baseCoordinate[i].pt2, PointCenter, N, distance, baseCoordinate[i].direction, PicCenter, &tempX, &tempY, &tempPsize);
		
		if(success)
		{
			//��ƽ����߾���
			CString tmp;
			tmp.Format("x:%f��y:%f��PixelSize:%f\n", tempX, tempY, tempPsize);
			Outcome+=tmp;
			//printf("x:%f��y:%f��PixelSize:%f\n", tempX, tempY, tempPsize);
			CenterX += tempX;
			CenterY += tempY;
			Psize += tempPsize;
			SuccessUnitNum++;
			
			//����ʱ������ȹ������������������߼���ٶ�
			//if(SuccessUnitNum==5)
			//	break;
		}		
	}
	if(SuccessUnitNum==0)//û�м���һ�����õĵ�Ԫ
	{
		free(PointCenter);
		free(baseCoordinate);
		free(pDataGaussian);
		free(pData);
		free(typeImg);
		return 0;
	}
	
	center->x = (DWORD)(CenterX/SuccessUnitNum/0.01+0.5);
	center->y = (DWORD)(CenterY/SuccessUnitNum/0.01+0.5);
	(*pixelsize) = (int)(Psize/SuccessUnitNum/0.01+0.5);


	free(PointCenter);
	free(baseCoordinate);
	free(pDataGaussian);
	free(pData);
	free(typeImg);
	return 1;
}

#define S 23
#define T (0.15f)

void adaptiveThreshold_C(BYTE* input, int IMAGE_WIDTH, int IMAGE_HEIGHT, int IMAGE_WIDESTEP, BYTE* bin)
{
	unsigned long* integralImg = 0;
	int i, j;
	long sum=0;
	int count=0;
	int index;
	int x1, y1, x2, y2;
	int s2 = S/2;

	// create the integral image
	integralImg = (unsigned long*)malloc(IMAGE_WIDESTEP*IMAGE_HEIGHT*sizeof(unsigned long*));
	memset(integralImg,0,sizeof(unsigned long*)*IMAGE_WIDESTEP*IMAGE_HEIGHT);

	for (i=0; i<IMAGE_WIDTH; i++)
	{
		// reset this column sum
		sum = 0;

		for (j=0; j<IMAGE_HEIGHT; j++)
		{
			index = j*IMAGE_WIDESTEP+i;

			sum += input[index];
			if (i==0)
				integralImg[index] = sum;
			else
				integralImg[index] = integralImg[index-1] + sum;
		}
	}

	// perform thresholding
	for (i=0; i<IMAGE_WIDTH; i++)
	{
		for (j=0; j<IMAGE_HEIGHT; j++)
		{
			index = j*IMAGE_WIDESTEP+i;

			// set the SxS region
			x1=i-s2; x2=i+s2;
			y1=j-s2; y2=j+s2;

			// check the border
			if (x1 < 0) x1 = 0;
			if (x2 >= IMAGE_WIDTH) x2 = IMAGE_WIDTH-1;
			if (y1 < 0) y1 = 0;
			if (y2 >= IMAGE_HEIGHT) y2 = IMAGE_HEIGHT-1;

			count = (x2-x1)*(y2-y1);

			// I(x,y)=s(x2,y2)-s(x1,y2)-s(x2,y1)+s(x1,x1)
			sum = integralImg[y2*IMAGE_WIDESTEP+x2] -
				integralImg[y1*IMAGE_WIDESTEP+x2] -
				integralImg[y2*IMAGE_WIDESTEP+x1] +
				integralImg[y1*IMAGE_WIDESTEP+x1];

			if ((long)(input[index]*count) < (long)(sum*(1.0-T)))
				bin[index] = 0;
			else
				bin[index] = 255;
		}
	}

	free (integralImg);
}

void gaussianFilter(unsigned char* corrupted, unsigned char* smooth, int width, int height)
{
	int i,j,sum,index,m,n;
	int pos;
	int templates[25] = { 1, 4, 7, 4, 1, 
						  4, 16, 26, 16, 4, 
						  7, 26, 41, 26, 7,
						  4, 16, 26, 16, 4, 
						  1, 4, 7, 4, 1 };//273

	//int templates[49] = { 1, 4, 7, 10, 7, 4, 1,//34
	//					  4, 12, 26, 33, 26, 12, 4,//117
	//					  7, 26, 55, 71, 55, 26, 7,//247
	//					  10, 33, 71, 91, 71, 33, 10,//319
	//					  7, 26, 55, 71, 55, 26, 7,//247
	//					  4, 12, 26, 33, 26, 12, 4,//117
	//					  1, 4, 7, 10, 7, 4, 1};//34;1115
	
	memcpy ( smooth, corrupted, width*height*sizeof(unsigned char) );
	for (j=0;j<height;j++)
	{
		for (i=0;i<width;i++)
		{
			sum = 0;
			index = 0;
			for ( m=j-2; m<j+3; m++)
			{
				for (n=i-2; n<i+3; n++)
				{
					pos = m*width + n;
					if(pos<0)
						pos=0;
					else if( pos>(height-1)*(width-1) )
						pos = (height-1)*(width-1);
					sum += corrupted[pos] * templates[index++];
				}
			}
			sum /= 273;
			if (sum > 255)
				sum = 255;
			smooth [ j*width+i ] = sum;
		}
	}
}