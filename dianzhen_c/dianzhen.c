#include "dianzhen.h"
#include "region.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <memory.h>

#include <time.h>
long t_start,t_end;

int g_ImgX = 0;
int g_ImgY = 0;
double g_distance = 0.0;//去除掉边界
bool IsPointInImg(CPoint pt)
{
	if(pt.x>=g_distance && pt.y>=g_distance && pt.x<g_ImgX-g_distance && pt.y<g_ImgY-g_distance)
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
	int i = 0;
	if(!IsPointInImg(pt))//如果该点不在图中
		return -2;
	for(i = 0; i < N; i++)
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

	double sigma = 1.8;//适应度参数，该参数应该小于2
	//const int des5_1Num = 5;//找5个1
	int pt5Idx[des5_1Num]={0};
	int i,j,ptNum,tep,ab_i;
	int curridx,nextInResion;
	int delta_x,delta_y;
	int zeroNum;
	
	double tempDis;
	CPoint beforePt,nextPt,pt1,pt2;
	CPoint a[ab_num],b[ab_num];
	CPoint tmpPta,tmpPtb;
	
	bool bFindBaseCoordinate;
	int BaseCoordinate_i;
	int a0_Region,b0_Region;
	
	for(i = 0; i < N; i++)
	{
		for(j = i+1; j < N; j++)
		{
			//if(i == 22)
			//	int mydebug = 1;
			//int endpt1_idx,endpt2_idx;//找到5个连续1的点中的两个端点
			//int curridx,nextInResion;
			tempDis = Get2PointDis(center[i], center[j]);
			if(tempDis > sigma*distance)//此点周围没有连续点，舍弃
				continue;			

			endpt1_idx = i;
			endpt2_idx = j;

			delta_x = center[j].x - center[i].x;
			delta_y = center[j].y - center[i].y;
			pt5Idx[0] = i;
			pt5Idx[1] = j;
			ptNum = 0;

			//先找左侧的目标点1
			//CPoint beforePt;
			beforePt.x = center[i].x - delta_x;
			beforePt.y = center[i].y - delta_y;
			curridx = i;
			nextInResion = PointInRegion(center, N, beforePt,distance);
			while(nextInResion>=0 && ptNum < 3)//下一坐标的点是目标点1
			{
				endpt1_idx = nextInResion;
				pt5Idx[ptNum+2] = nextInResion;
				ptNum++;
				beforePt.x = 2*center[nextInResion].x - center[curridx].x;
				beforePt.y = 2*center[nextInResion].y - center[curridx].y;

				curridx = nextInResion;
				nextInResion = PointInRegion(center, N, beforePt,distance);
			}

			//后找右侧的目标点1
			//CPoint nextPt;
			nextPt.x = center[j].x+delta_x;
			nextPt.y = center[j].y+delta_y;
			curridx = j;
			nextInResion = PointInRegion(center, N, nextPt,distance);
			while(nextInResion>=0 && ptNum < des5_1Num-2)//下一坐标的点是目标点1
			{
				endpt2_idx = nextInResion;
				pt5Idx[ptNum+2] = nextInResion;
				ptNum++;
				nextPt.x = 2*center[nextInResion].x - center[curridx].x;
				nextPt.y = 2*center[nextInResion].y - center[curridx].y;

				curridx = nextInResion;
				nextInResion = PointInRegion(center, N, nextPt,distance);
			}

			if(ptNum!=des5_1Num-2)//必须满足是5个连续的1才可能是坐标
				continue;

			//判断端点是否已经被找出过
			bFindBaseCoordinate = false;
			for(BaseCoordinate_i = 0; BaseCoordinate_i < BaseCoordinateNum; BaseCoordinate_i++)
			{
				if(baseCoordinate[BaseCoordinate_i].pt1 == endpt1_idx || baseCoordinate[BaseCoordinate_i].pt1 == endpt2_idx)
				{
					bFindBaseCoordinate = true;
					break;
				}
			}
			if(bFindBaseCoordinate)
				continue;


			//判断两个端点两侧是否为0，若不是则这5个点不是坐标
			delta_4x = (center[endpt2_idx].x - center[endpt1_idx].x)/4.0;
			delta_4y = (center[endpt2_idx].y - center[endpt1_idx].y)/4.0;
			//CPoint pt1;
			pt1.x = (int)(center[endpt1_idx].x-delta_4x);
			pt1.y = (int)(center[endpt1_idx].y-delta_4y);
			if(PointInRegion(center, N, pt1, distance) != -1)//判断是否为0，若不是则这5个点不是坐标,舍弃
				continue;
			//CPoint pt2;
			pt2.x = (int)(center[endpt2_idx].x+delta_4x);
			pt2.y = (int)(center[endpt2_idx].y+delta_4y);
			tep = PointInRegion(center, N, pt2, distance);
			if(PointInRegion(center, N, pt2, distance) != -1)//判断是否为0，若不是则这5个点不是坐标,舍弃
				continue;

			//在5个点两侧pt1,pt2是0的基础上，判断pt1,pt2上下是否为10/11
			//const int ab_num = 3;
			//CPoint a[ab_num],b[ab_num];//(above,below)
			for(ab_i = 0; ab_i < ab_num; ab_i++)
			{
				a[ab_i].x = (int)(pt1.x + (ab_i+1)*delta_4y);
				a[ab_i].y = (int)(pt1.y - (ab_i+1)*delta_4x);

				b[ab_i].x = (int)(pt2.x - (ab_i+1)*delta_4y);
				b[ab_i].y = (int)(pt2.y + (ab_i+1)*delta_4x);
			}

			if(PointInRegion(center, N, a[1], distance)<0 || PointInRegion(center, N, b[1], distance)<0)//点a1,b1必须都为1，否则舍弃
				continue;
			if(PointInRegion(center, N, a[2], distance)>=0 || PointInRegion(center, N, b[2], distance)>=0)//点a2,b2必须都为0，否则舍弃
				continue;


			//判断是否有a2,b2的右边是否有5个0
			zeroNum = 1;
			//CPoint tmpPta,tmpPtb;

			while(zeroNum<=6)
			{
				tmpPta.x = a[2].x+zeroNum*delta_4x;
				tmpPta.y = a[2].y+zeroNum*delta_4y;

				tmpPtb.x = b[2].x+zeroNum*delta_4x;
				tmpPtb.y = b[2].y+zeroNum*delta_4y;
				if(PointInRegion(center, N, tmpPta, distance)>=0 || PointInRegion(center, N, tmpPtb, distance)>=0)//如果该位置不是0，跳出
					break;
				zeroNum++;

			}
			if(zeroNum<6)
				continue;


			a0_Region = PointInRegion(center, N, a[0], distance);
			b0_Region = PointInRegion(center, N, b[0], distance);

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
			if(BaseCoordinateNum == 5)
				return BaseCoordinateNum;
		}
	}

	return BaseCoordinateNum;
}

//求向量的内积
double DotProduct(CPoint a, CPoint b)
{
	return a.x*b.x + a.y*b.y;
}


bool GetXyCoordinate(int endPt1, int endPt2, CPoint* PointCenter, int N, double distance, bool direction, CPoint PicCenter, double* tempX, double* tempY, double* PSize)
{
	CPoint BasePt = PointCenter[endPt1];
	double delta_x = (PointCenter[endPt2].x - PointCenter[endPt1].x)/4.0;
	double delta_y = (PointCenter[endPt2].y - PointCenter[endPt1].y)/4.0;

	//const int xynum = 10;
	CPoint xPt[xynum],yPt[xynum];	
	CPoint StartPt,EndPt;
	CPoint V_X,V_Y;
	CPoint UnitCenter,V_U2P;
	
	double xyDotP,V_Xmode,V_Ymode,x_shadow,y_shadow,PixelsDis_2Point;

	int x_ratio,y_ratio;//两层系数不同
	int xOutAarry[xynum],yOutAarry[xynum];
	int xytmp,xOut,yOut;
	
	
	
	int i;
	for(i = 0; i < xynum; i++)
	{	
		x_ratio = i>4?1:2;
		xPt[i].x = (int)(BasePt.x + x_ratio*delta_y + (i%5)*delta_x);
		xPt[i].y = (int)(BasePt.y - x_ratio*delta_x + (i%5)*delta_y);

		y_ratio = i>4?2:1;
		yPt[i].x = (int)(BasePt.x - y_ratio*delta_y + (i%5)*delta_x);
		yPt[i].y = (int)(BasePt.y + y_ratio*delta_x + (i%5)*delta_y);

	}

	//如果四个端点不在图像中，则跳出
	if(!IsPointInImg(xPt[0]) || !IsPointInImg(xPt[4]) || !IsPointInImg(yPt[5]) || !IsPointInImg(yPt[9]))
		return false;

	//int xOutAarry[xynum],yOutAarry[xynum];
	for(i = 0; i < xynum; i++)
	{
		xOutAarry[i] = PointInRegion(PointCenter, N, xPt[i], distance)>=0?1:0;
		yOutAarry[i] = PointInRegion(PointCenter, N, yPt[i], distance)>=0?1:0;
	}

	if(!direction)//如果方向不对，必须调换方向
	{
		xytmp;
		for(i = 0; i < xynum; i++)
		{
			xytmp = xOutAarry[i];
			xOutAarry[i] = yOutAarry[9-i];
			yOutAarry[9-i] = xytmp; 
		}
	}

	xOut = 0;
	yOut = 0;

	for(i = 0; i < xynum; i++)
	{
		xOut = (xOut<<1) | xOutAarry[i];
		yOut = (yOut<<1) | yOutAarry[i];
	}
	//printf("坐标x:%d,y:%d\t\t",xOut, yOut );

	//计算原始坐标系x，y坐标轴方向
	//CPoint StartPt,EndPt;

	StartPt.x = direction?xPt[7].x:yPt[2].x;
	StartPt.y = direction?xPt[7].y:yPt[2].y;
	EndPt.x = direction?xPt[9].x:yPt[0].x;
	EndPt.y = direction?xPt[9].y:yPt[0].y;

	//CPoint V_X;//={(), (EndPt.y-StartPt.y)};//原始x坐标向量
	V_X.x=EndPt.x-StartPt.x;
	V_X.y=EndPt.y-StartPt.y;

	StartPt.x = direction?xPt[7].x:yPt[2].x;
	StartPt.y = direction?xPt[7].y:yPt[2].y;
	EndPt.x = direction?yPt[2].x:xPt[7].x;
	EndPt.y = direction?yPt[2].y:xPt[7].y;


	//CPoint V_Y;//={(EndPt.x-StartPt.x), (EndPt.y-StartPt.y)};//原始y坐标向量
	V_Y.x=EndPt.x-StartPt.x;
	V_Y.y=EndPt.y-StartPt.y;

	xyDotP = DotProduct(V_X, V_Y);//求xy坐标向量的内积，因为x,y坐标垂直，期望值应该为0，大多数能符合说明方向是正确的

	//求xy坐标向量的模
	V_Xmode = sqrt(1.0*(V_X.x*V_X.x + V_X.y*V_X.y));
	V_Ymode = sqrt(1.0*(V_Y.x*V_Y.x + V_Y.y*V_Y.y));


	//检测单元的中心点
	//CPoint UnitCenter;
	UnitCenter.x = (xPt[7].x + yPt[2].x)/2;
	UnitCenter.y = (xPt[7].y + yPt[2].y)/2;


	//CPoint V_U2P;//检测单元的中心点与图片中心组成的向量，由检测单元的中心指向图片中心
	V_U2P.x = PicCenter.x - UnitCenter.x;
	V_U2P.y = PicCenter.y - UnitCenter.y;



	//计算向量V_U2P在x，y方向的投影
	x_shadow = DotProduct(V_X, V_U2P)/V_Xmode;
	y_shadow = DotProduct(V_Y, V_U2P)/V_Ymode;

	PixelsDis_2Point = sqrt(delta_x*delta_x + delta_y*delta_y);
	(*PSize) = REALDIS/PixelsDis_2Point ;//像素代表的实际距离

	(*tempX) = (*PSize) * x_shadow + xOut*1.98;//x方向的偏移量
	(*tempY) = (*PSize) * y_shadow + yOut*1.98;//y方向的偏移量

	return true;
}

int small_gaussian[4][4] =
{
	{1024}, 
	{512, 256}, 
	{384 , 256, 64},
	{288, 224, 112, 32}
};


//多通道图像高斯平滑
void MyGaussian(unsigned char *pSrcData, unsigned char *pDesData, int Width, int Height)
{
	int i,j,k,m;
	//int Width=pSrc->width;
	//int Height=pSrc->height;
	int nCh=1;
	int widestep=Width;
	unsigned char * Data=pSrcData;	
	unsigned char * sData;
	int SmoothRadio=3;
	int pos,rpos;
	int Temp;
	int *  pGauTepData=NULL;
	unsigned char * rfilter;

	//根据SMOOTHRADIO判断是哪个数组
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

	rfilter=(unsigned char*)malloc(sizeof(unsigned char)*Width*Height);
	sData=rfilter;
	//行滤波
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

void BilinearRowFilter(unsigned char* src, long* dst, int len, long* leftIdx, long* rightIdx, long* weight, int shift ,int nch) 
{ 
	int i,j;	
	for(i=0;i<len;i++)
	{
		for (j=0;j<nch;j++)
		{
			*dst++	=((1<<shift) - weight[i])*src[leftIdx[i]*nch+j] + weight[i]*src[rightIdx[i]*nch+j]; 
		}
	}
}

void BiLinearInsert(unsigned char * pSrc, int sWidth, int sHeight, unsigned char * pDes, int dWidth, int dHeight ,int nch)
{
	int linewidth = sWidth;//(sWidth*nch+3)/4*4;//通道数为1，不考虑widthStep优化，即传进来的数据全部是有效数据
	int dlinewidth = dWidth;//(dWidth*nch+3)/4*4;//通道数为1，不考虑widthStep优化，即传进来的数据全部是有效数据

	float dx=1.0f*sWidth/dWidth;
	float dy=1.0f*sHeight/dHeight;
	int i,j,m;
	int shift=10;
	
	long *rowBuf1,*rowBuf2,*horWeight,*horLeftIdx,*horRightIdx;
	int preVerDownIdx,preVerUpIdx;
	long *upLinePtr,*downLinePtr,*tempPtr;
	
	float pos;
	int verUpIdx,verDownIdx;
	long verWeight;
	unsigned char* _ptr;

	rowBuf1=(long*)malloc(sizeof(long)*dlinewidth);//new long[dlinewidth];
	rowBuf2=(long*)malloc(sizeof(long)*dlinewidth);;

	horWeight = (long*)malloc(sizeof(long)*dWidth);//new long[dWidth]; 
	horLeftIdx = (long*)malloc(sizeof(long)*dWidth);//new long[dWidth]; 
	horRightIdx = (long*)malloc(sizeof(long)*dWidth);//new long[dWidth]; 
	
	preVerUpIdx = -1; 
	preVerDownIdx = -1; 
	upLinePtr = rowBuf1; 
	downLinePtr = rowBuf2; 
	tempPtr;

	for(i = 0; i < dWidth; i++)
	{ 
		pos = (i + 0.5f) * dx; 
		horLeftIdx[i] = (int)(MAX(pos - 0.5f, 0)); 
		horRightIdx[i] = (int)(MIN(pos + 0.5f, sWidth-1)); 
		horWeight[i] = (long) (fabs(pos - 0.5f - horLeftIdx[i]) * (1<<shift)); 
	} 


	for(j = 0; j < dHeight; j++)
	{ 
		pos = (j + 0.5f) * dy; 
		verUpIdx = (int)(MAX(pos - 0.5f, 0));  
		verDownIdx = (int)(MIN(pos + 0.5f, sHeight-1)); 
		verWeight = (long) (fabs(pos - 0.5f - verUpIdx) * (1<<shift)); 

		if(verUpIdx == preVerUpIdx && verDownIdx == preVerDownIdx)
		{ 
			;
		} 
		else if(verUpIdx == preVerDownIdx)  
		{ 
			
			tempPtr=upLinePtr;
			upLinePtr=downLinePtr;
			downLinePtr=tempPtr;
			BilinearRowFilter(pSrc + linewidth*verDownIdx, downLinePtr, dWidth, horLeftIdx, horRightIdx, horWeight, shift,nch); 
		}
		else
		{ 
			
			BilinearRowFilter(pSrc + linewidth*verUpIdx,   upLinePtr, dWidth, horLeftIdx, horRightIdx, horWeight, shift,nch);
			
			BilinearRowFilter(pSrc + linewidth*verDownIdx, downLinePtr, dWidth, horLeftIdx, horRightIdx, horWeight, shift,nch); 
		} 
		
		_ptr = pDes + j*dlinewidth; 
		for(i = 0; i < dWidth; i++)
		{ 
			for (m=0;m<nch;m++)
			{
				*_ptr++=(unsigned char) ( (((1<<shift) - verWeight)*upLinePtr[i*nch+m] + verWeight*downLinePtr[i*nch+m]) >> (2*shift) ); 
			}
		} 

		preVerUpIdx = verUpIdx; 
		preVerDownIdx = verDownIdx; 
	} 

	free(rowBuf1);
	free(rowBuf2);
	free(horWeight);
	free(horLeftIdx);
	free(horRightIdx);
}

int process(unsigned char * y2ybuf, int width, int height, Axis *center, int *pixelsize)
{

	unsigned char *pDataGaussian,*resizebuf;
	int * typeImg;
	int i,j,N,True_region_num,w,h;
	extern Region region[1000];
	Region* TrueRegion;
	CPoint* PointCenter;
	double distance,boxRatio,tempDis;
	int topBox,bottomBox,leftBox,rightBox;
	CPoint PicCenter;
	BaseCoordinate* baseCoordinate;
	int BaseCoordinateNum,SuccessUnitNum;
	double CenterX,CenterY,Psize,tempX,tempY,tempPsize;
	bool success;
	
	//0.图像缩放
	double PicScale;
	int resize_width,resize_height;

	PicScale = 300.0/width;
	
	resize_width = 300;
	resize_height = (int)(PicScale*height);
	t_start = clock();
	resizebuf = (unsigned char*)malloc(sizeof(unsigned char)*(resize_width*resize_height));
	BiLinearInsert(y2ybuf, width, height, resizebuf, resize_width, resize_height, 1);
	t_end = clock();
	printf("BiLinearInsert:%d\n",t_end-t_start);
	t_start = clock();
	g_ImgX = resize_width;
	g_ImgY = resize_height;

	//1.高斯滤波
	pDataGaussian = (unsigned char*)malloc(sizeof(unsigned char)*(resize_width*resize_height));

	memcpy(pDataGaussian, resizebuf, sizeof(unsigned char)*(resize_width*resize_height));

	MyGaussian(resizebuf, pDataGaussian, resize_width, resize_height);
	t_end = clock();
	printf("MyGaussian:%d\n",t_end-t_start);
	//2.二值化
	t_start=clock();
	adaptiveThreshold_C(pDataGaussian, resize_width, resize_height, resize_width, resizebuf);
	t_end = clock();
	printf("adaptiveThreshold_C:%d\n",t_end-t_start);
	//3.寻找连通区域
	t_start=clock();
	typeImg=(int*)malloc(sizeof(int)*resize_height*resize_width);
	for(i=0;i<resize_height;i++)
	{
		for(j=0;j<resize_width;j++)
		{
			typeImg[i*resize_width+j]=-1;
		}
	}
	N=KeyGrow(resizebuf,resize_width,resize_height,typeImg);	

	

	True_region_num = N;
	for (i=0;i<N;i++)
	{
		//宽度比较大的，并且占空比比较大
		w=region[i].right-region[i].left+1;
		h=region[i].bottom-region[i].top+1;
		//if ( region[i].PointNum<=5)
		//{
		//	region[i].IsOK=false;
		//	True_region_num--;
		//}

		//double tmp_ratio = w*1.0/h;
		//if(tmp_ratio>2||tmp_ratio<0.5)
		//{
		//	region[i].IsOK=false;
		//	True_region_num--;
		//}

		//太高的太长的并且占空比比较大的去掉
		//if (w<=2&&h<=2)
		//{
		//	region[i].IsOK=false;
		//	True_region_num--;
		//}
	}
	TrueRegion = (Region*)malloc(sizeof(Region)*N);
	True_region_num = 0;
	for(i = 0; i< N; i++)
	{
		if(region[i].IsOK)
		{
			TrueRegion[True_region_num].bottom = region[i].bottom;
			TrueRegion[True_region_num].left = region[i].left;
			TrueRegion[True_region_num].right = region[i].right;
			TrueRegion[True_region_num].top = region[i].top;

			True_region_num++;
		}
	}

	N = True_region_num;

	//4.通过找最小的两点距离确定栅格距离distance

	PointCenter = (CPoint*)malloc(sizeof(CPoint)*N);//求连通域的最小外接矩阵的中心
	for(i = 0; i < N; i++)
	{
		PointCenter[i].x = (TrueRegion[i].left+TrueRegion[i].right)>>1;
		PointCenter[i].y = (TrueRegion[i].top+TrueRegion[i].bottom)>>1;

	}


	distance = 0xFFFFFFF;
	boxRatio = 0.2;//找最短距离时去除在图片四周的点，防止找到错误的diastance
	leftBox = boxRatio * resize_width;
	rightBox = resize_width - leftBox;
	topBox = boxRatio * resize_height;
	bottomBox = resize_height - topBox;
	for(i = 0; i < N; i++)
	{
		if(PointCenter[i].x < leftBox || PointCenter[i].x > rightBox || PointCenter[i].y < topBox || PointCenter[i].y > bottomBox)
			continue;
		for(j = i+1; j < N; j++)
		{
			if(PointCenter[j].x < leftBox || PointCenter[j].x > rightBox || PointCenter[j].y < topBox || PointCenter[j].y > bottomBox)
				continue;
			tempDis = Get2PointDis(PointCenter[i], PointCenter[j]);
			if(tempDis < distance)
			{
				distance = tempDis;
			}
		}
	}

	g_distance = sqrt(distance);//去除掉边界
	//CPoint PicCenter;//={width/2, height/2};
	PicCenter.x=resize_width/2;
	PicCenter.y=resize_height/2;
	t_end = clock();
	printf("findDis:%d\n",t_end-t_start);
	//5寻找基坐标
	t_start=clock();
	baseCoordinate = (BaseCoordinate*)malloc(sizeof(BaseCoordinate)*N/8);//每个栅格最少有8个1
	BaseCoordinateNum = FindBaseCoordinate(PointCenter, N, distance, baseCoordinate);
	t_end=clock();
	printf("FindBaseCoordinate:%d\n",t_end-t_start);
	SuccessUnitNum = 0;

	center->x = 0;
	center->y = 0;
	CenterX = 0;
	CenterY = 0;
	Psize = 0;
	tempX,tempY,tempPsize;
	t_start=clock();
	for(i = 0; i < BaseCoordinateNum; i++)
	{	
		
		success = GetXyCoordinate(baseCoordinate[i].pt1, baseCoordinate[i].pt2, PointCenter, N, distance, baseCoordinate[i].direction, PicCenter, &tempX, &tempY, &tempPsize);

		if(success)
		{
			//求平均提高精度
			CenterX += tempX;
			CenterY += tempY;
			Psize += tempPsize;
			SuccessUnitNum++;

			////测试时如果精度够，可以提高跳出以提高检测速度
			//if(SuccessUnitNum==5)
			//	break;
		}		
	}
	t_end=clock();
	printf("GetXyCoordinate:%d\n",t_end-t_start);

	if(SuccessUnitNum==0)//没有检测出一个可用的单元
	{
		free(PointCenter);
		free(baseCoordinate);
		free(pDataGaussian);
		free(resizebuf);
		free(typeImg);
		return 0;
	}

	center->x = (DWORD)(CenterX/SuccessUnitNum/0.01+0.5);
	center->y = (DWORD)(CenterY/SuccessUnitNum/0.01+0.5);
	(*pixelsize) = (int)(Psize/SuccessUnitNum/0.001*PicScale+0.5);


	free(PointCenter);
	free(baseCoordinate);
	free(pDataGaussian);
	free(resizebuf);
	free(typeImg);
	free(TrueRegion);

	return 1;
}

#define S 23
#define T (0.15f)

void adaptiveThreshold_C(unsigned char* input, int IMAGE_WIDTH, int IMAGE_HEIGHT, int IMAGE_WIDESTEP, unsigned char* bin)
{
	
	unsigned long* integralImg = 0;
	int i, j;
	long sum=0;
	int count=0;
	int index;
	int x1, y1, x2, y2;
	int s2 = S/2;
	t_start = clock();
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
	/*t_end = clock();
	printf("integralImg:%d\n", t_end-t_start);
	t_start = clock();*/
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
	t_end = clock();
	//printf("Threshold:%d\n", t_end-t_start);
}
