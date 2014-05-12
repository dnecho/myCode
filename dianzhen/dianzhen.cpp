#include "stdafx.h"
int g_ImgX = 0;
int g_ImgY = 0;
inline bool IsPointInImg(CPoint pt)
{
	if(pt.x>=0 && pt.y>=0 && pt.x<g_ImgX && pt.y<g_ImgY)
		return true;
	return false;
}

inline double Get2PointDis(CPoint a, CPoint b)
{
	return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
}

int PointInRegion(CPoint* center, int N, CPoint pt, double dis)
{
	double d = dis*0.5;
	if(!IsPointInImg(pt))//如果该点不在图中
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

	double sigma = 1.95;//适应度参数，该参数应该小于2
	const int des5_1Num = 5;//找5个1
	int pt5Idx[des5_1Num]={0};
	for(int i = 0; i < N; i++)
	{
		for(int j = i+1; j < N; j++)
		{
			//if(i == 22)
			//	int mydebug = 1;
			//int endpt1_idx,endpt2_idx;//找到5个连续1的点中的两个端点
			int curridx,nextInResion;
			double tempDis = Get2PointDis(center[i], center[j]);
			if(tempDis > sigma*distance)//此点周围没有连续点，舍弃
				continue;			

			endpt1_idx = i;
			endpt2_idx = j;

			int delta_x = center[j].x - center[i].x;
			int delta_y = center[j].y - center[i].y;
			pt5Idx[0] = i;
			pt5Idx[1] = j;
			int ptNum = 0;

			//先找左侧的目标点1
			CPoint beforePt;
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
			CPoint nextPt;
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


			//判断两个端点两侧是否为0，若不是则这5个点不是坐标
			delta_4x = (center[endpt2_idx].x - center[endpt1_idx].x)/4.0;
			delta_4y = (center[endpt2_idx].y - center[endpt1_idx].y)/4.0;
			CPoint pt1;
			pt1.x = (int)(center[endpt1_idx].x-delta_4x);
			pt1.y = (int)(center[endpt1_idx].y-delta_4y);
			if(PointInRegion(center, N, pt1, distance) != -1)//判断是否为0，若不是则这5个点不是坐标,舍弃
				continue;
			CPoint pt2;
			pt2.x = (int)(center[endpt2_idx].x+delta_4x);
			pt2.y = (int)(center[endpt2_idx].y+delta_4y);
			int tep = PointInRegion(center, N, pt2, distance);
			if(PointInRegion(center, N, pt2, distance) != -1)//判断是否为0，若不是则这5个点不是坐标,舍弃
				continue;

			//在5个点两侧pt1,pt2是0的基础上，判断pt1,pt2上下是否为10/11
			const int ab_num = 2;
			CPoint a[ab_num],b[ab_num];//(above,below)
			for(int ab_i = 0; ab_i < ab_num; ab_i++)
			{
				a[ab_i].x = (int)(pt1.x + delta_4y);
				a[ab_i].y = (int)(pt1.y - delta_4x);

				b[ab_i].x = (int)(pt1.x - delta_4y);
				b[ab_i].y = (int)(pt1.y + delta_4x);
			}

			if(PointInRegion(center, N, a[1], distance)<0 && PointInRegion(center, N, b[1], distance)<0)//点a1,b1必须都为1，否则舍弃
				continue;
			int a0_Region = PointInRegion(center, N, a[0], distance);
			int b0_Region = PointInRegion(center, N, b[0], distance);
			//char txt[255];
			//extern Region region[1000];
			if(a0_Region==-1 && b0_Region>=0 )
			{
				//for(int pt_i = 0; pt_i < des5_1Num; pt_i++)
				//{
				//	sprintf_s(txt, "%d", pt_i);

				//	cvRectangle(g_pImg, cvPoint(region[pt5Idx[pt_i]].left, region[pt5Idx[pt_i]].top), cvPoint(region[pt5Idx[pt_i]].right, region[pt5Idx[pt_i]].bottom), cvScalar(255,0,0));
				//	//cout<<pt5Idx[pt_i]<<"\t";
				//}
				//cout<<"BaseCoordinateNum:"<<BaseCoordinateNum<<"\tendpt1:"<<endpt1_idx<<"\t"<<"endpt2:"<<endpt2_idx<<endl;
				baseCoordinate[BaseCoordinateNum].pt1 = endpt1_idx;
				baseCoordinate[BaseCoordinateNum].pt2 = endpt2_idx;
				baseCoordinate[BaseCoordinateNum].direction = true;
				BaseCoordinateNum++;
				//direction = true;
				//return true;
			}
			else if(a0_Region>0 && b0_Region==-1 )
			{
				//for(int pt_i = 0; pt_i < des5_1Num; pt_i++)
				//{
				//	sprintf_s(txt, "%d", pt_i);

				//	cvRectangle(g_pImg, cvPoint(region[pt5Idx[pt_i]].left, region[pt5Idx[pt_i]].top), cvPoint(region[pt5Idx[pt_i]].right, region[pt5Idx[pt_i]].bottom), cvScalar(255,0,0));
				//	//cout<<pt5Idx[pt_i]<<"\t";					
				//}
				//cout<<"BaseCoordinateNum:"<<BaseCoordinateNum<<"\tendpt1:"<<endpt1_idx<<"\t"<<"endpt2:"<<endpt2_idx<<endl;
				baseCoordinate[BaseCoordinateNum].pt1 = endpt1_idx;
				baseCoordinate[BaseCoordinateNum].pt2 = endpt2_idx;
				baseCoordinate[BaseCoordinateNum].direction = false;
				BaseCoordinateNum++;
				//direction = false;
				//return true;
			}				
		}
	}

	return BaseCoordinateNum;
}

bool GetXyCoordinate(IplImage* pOutImg, int endPt1, int endPt2, CPoint* center, int N, double distance, bool direction)
{
	cout<<"1:"<<endPt1<<"\t2:"<<endPt2<<endl;
	CPoint BasePt = center[endPt1];
	double delta_x = (center[endPt2].x - center[endPt1].x)/4.0;
	double delta_y = (center[endPt2].y - center[endPt1].y)/4.0;

	const int xynum = 10;
	CPoint xPt[xynum],yPt[xynum];
	int x_ratio,y_ratio;//两层系数不同
	for(int i = 0; i < xynum; i++)
	{	
		x_ratio = i>4?1:2;
		xPt[i].x = (int)(BasePt.x + x_ratio*delta_y + (i%5)*delta_x);
		xPt[i].y = (int)(BasePt.y - x_ratio*delta_x + (i%5)*delta_y);

		y_ratio = i>4?2:1;
		yPt[i].x = (int)(BasePt.x - y_ratio*delta_y + (i%5)*delta_x);
		yPt[i].y = (int)(BasePt.y + y_ratio*delta_x + (i%5)*delta_y);

		//sprintf_s(txt, "%d", i);
		//cvPutText(pSrc3C, txt, cvPoint(xPt[i].x, xPt[i].y), &font,cvScalar(0,255,255));
		//sprintf_s(txt, "%d", i);
		//cvPutText(pSrc3C, txt, cvPoint(yPt[i].x, yPt[i].y), &font,cvScalar(0,255,255));
	}

	//如果四个端点不在图像中，则跳出
	if(!IsPointInImg(xPt[0]) || !IsPointInImg(xPt[4]) || !IsPointInImg(yPt[5]) || !IsPointInImg(yPt[9]))
		return false;

	cvLine(pOutImg, cvPoint(xPt[0].x,xPt[0].y), cvPoint(xPt[4].x,xPt[4].y), cvScalar(0,255,0));
	cvLine(pOutImg, cvPoint(yPt[9].x,yPt[9].y), cvPoint(xPt[4].x,xPt[4].y), cvScalar(0,255,0));
	cvLine(pOutImg, cvPoint(yPt[9].x,yPt[9].y), cvPoint(yPt[5].x,yPt[5].y), cvScalar(0,255,0));
	cvLine(pOutImg, cvPoint(xPt[0].x,xPt[0].y), cvPoint(yPt[5].x,yPt[5].y), cvScalar(0,255,0));
	int xOutAarry[xynum],yOutAarry[xynum];
	for(int i = 0; i < xynum; i++)
	{
		xOutAarry[i] = PointInRegion(center, N, xPt[i], distance)>=0?1:0;
		yOutAarry[i] = PointInRegion(center, N, yPt[i], distance)>=0?1:0;
	}

	if(!direction)//如果方向不对，必须调换方向
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

	CvFont font; 
	cvInitFont(&font,CV_FONT_HERSHEY_PLAIN ,1,1,0,1);

	CvPoint xPos,yPos;

	xPos.x = direction?xPt[7].x:yPt[2].x;
	xPos.y = direction?xPt[7].y:yPt[2].y;
	yPos.x = direction?yPt[2].x:xPt[7].x;
	yPos.y = direction?yPt[2].y:xPt[7].y;

	drawArrow(pOutImg, xPos, yPos, 5, 30, cvScalar(0,255,0));

	char txt[255];
	sprintf_s(txt, "%s%d", "x:", xOut);
	cvPutText(pOutImg, txt, xPos, &font,cvScalar(0,0,255));
	sprintf_s(txt, "%s%d", "y:", yOut);
	cvPutText(pOutImg, txt, yPos, &font,cvScalar(0,0,255));

	return true;
}

void Detect(IplImage* pSrc, IplImage* pSrc3C)
{
	int width = pSrc->width;
	int height = pSrc->height;
	g_ImgX = width;
	g_ImgY = height;

	//1.高斯滤波
	cvSmooth(pSrc, pSrc, CV_GAUSSIAN, 7);

	//2.二值化
	IplImage* pThresholdImg = cvCloneImage(pSrc);
	adaptiveThreshold(pSrc, pThresholdImg);
	unsigned char* pData = (unsigned char*)pThresholdImg->imageData;

	//3.寻找连通区域
	int * typeImg=new int[height*width];
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			typeImg[i*width+j]=-1;
		}
	}
	int N=KeyGrow(pData,width,height,typeImg);

	CvFont font; 

	cvInitFont(&font,CV_FONT_HERSHEY_PLAIN ,0.5,0.5,0,1);  


	extern Region region[1000];
	char txt[255]="";
	//for(int i=0;i<N;i++)
	//{
	//	if (region[i].IsOK)
	//	{	
	//		sprintf_s(txt, "%d", i);
	//		
	//		cvRectangle(pSrc3C, cvPoint(region[i].left, region[i].top), cvPoint(region[i].right, region[i].bottom), cvScalar(255,0,0));
	//		//cvPutText(pSrc3C, txt, cvPoint(region[i].left, region[i].top), &font,cvScalar(0,255,255));  
	//	}
	//}
	//cvNamedWindow("pSrc3C", 0);
	//cvShowImage("pSrc3C", pSrc3C);
	//cvWaitKey();
	//4.通过找最小的两点距离确定栅格距离distance	
	CPoint* center = new CPoint[N];//求连通域的最小外接矩阵的中心
	for(int i = 0; i < N; i++)
	{
		center[i].x = (region[i].left+region[i].right)>>1;
		center[i].y = (region[i].top+region[i].bottom)>>1;
	}

	double distance = 0xFFFFFFF;
	for(int i = 0; i < N; i++)
	{
		for(int j = i+1; j < N; j++)
		{
			double tempDis = Get2PointDis(center[i], center[j]);
			if(tempDis < distance)
				distance = tempDis;

		}
	}


	//5寻找基坐标
	BaseCoordinate* baseCoordinate = new BaseCoordinate[N/8];//每个栅格最少有8个1
	int BaseCoordinateNum = FindBaseCoordinate(center, N, distance, baseCoordinate);	
	for(int i = 0; i < BaseCoordinateNum; i++)
	{
		GetXyCoordinate(pSrc3C, baseCoordinate[i].pt1, baseCoordinate[i].pt2, center, N, distance, baseCoordinate[i].direction);
	}
	
	delete []center;
	delete []baseCoordinate;
	//cvNamedWindow("Src", 0);
	//cvShowImage("Src", pThresholdImg);
}

void drawArrow(IplImage* img, CvPoint pStart, CvPoint pEnd, int len, int alpha,  CvScalar& color, int thickness, int lineType)
{
	  const double PI = 3.1415926;    
	  Point arrow;    
	  //计算 θ 角（最简单的一种情况在下面图示中已经展示，关键在于 atan2 函数，详情见下面）   
	  double angle = atan2((double)(pStart.y - pEnd.y), (double)(pStart.x - pEnd.x));
	  cvLine(img, pStart, pEnd, color, thickness, lineType);
	  //line(img, pStart, pEnd, color, thickness, lineType);   
	  //计算箭角边的另一端的端点位置（上面的还是下面的要看箭头的指向，也就是pStart和pEnd的位置） 
	  arrow.x = (int)(pEnd.x + len * cos(angle + PI * alpha / 180));     
	  arrow.y = (int)(pEnd.y + len * sin(angle + PI * alpha / 180));  
	  //line(img, pEnd, arrow, color, thickness, lineType);
	  cvLine(img, pEnd, arrow, color, thickness, lineType);   
	  arrow.x = (int)(pEnd.x + len * cos(angle - PI * alpha / 180));     
	  arrow.y = (int)(pEnd.y + len * sin(angle - PI * alpha / 180));    
	  cvLine(img, pEnd, arrow, color, thickness, lineType);
 }