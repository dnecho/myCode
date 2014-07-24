#include "stdafx.h"
#include <memory.h>
#include <malloc.h>

Region region[REGION_NUM];

int KeyGrow(unsigned char * p, int w, int h,int * typeImg)
{
	int LineSize;
	unsigned char DealPixel;
	int * imgcopy;
	CPoint keytype[MaxPointNum];
	int N,nLEFT,nRight,nTOP,nBOTTOM;
	int i,j,n,typeNum;
	int x,y;
	
	if(p==0)
	{
		return 0;
	}

	//根据具体情况判断是否需要进行四字节对齐
	LineSize=(w+3)/4*4;
	DealPixel=0;

	imgcopy=(int*)malloc(sizeof(int)*LineSize*h);
	memset(imgcopy,0,sizeof(int)*LineSize*h);

	//const int MaxPointNum=20000;
	//CPoint keytype[MaxPointNum];

	//类别的小标
	N=-1;
	//int nLEFT,nRight,nTOP,nBOTTOM;
	for (j=0;j<h;j++)
	{
		for (i=0;i<w;i++)
		{
			//数组下标
			n=-1;
			//当前白点的数目
			typeNum=0;
			//如果当前的是白点，并且还没有处理过，将其作为一个种子点，向四周扩散
			if(imgcopy[j*LineSize+i]==0 && p[j*LineSize+i]==DealPixel )
			{
				typeImg[j*w+i]=N+1;

				n++;
				keytype[n].x=i;
				keytype[n].y=j;
				typeNum++;

				imgcopy[j*LineSize+i]=1;		
			}
			else
				continue;

			nLEFT=LineSize;
			nRight=0;
			nTOP=h;
			nBOTTOM=0;
			while (n>=0 && n<MaxPointNum)
			{
				//数组末尾的值相当于POP出来
				x=keytype[n].x;
				y=keytype[n].y;

				nLEFT=nLEFT<x?nLEFT:x;
				nRight=nRight>x?nRight:x;
				nTOP=nTOP<y?nTOP:y;
				nBOTTOM=nBOTTOM>y?nBOTTOM:y;

				//下标减一
				n--;

				//zuo
				if(x-1>=0&&p[y*LineSize+x-1]==DealPixel&&imgcopy[y*LineSize+x-1]==0)
				{
					typeImg[y*w+x-1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x-1;
					keytype[n].y=y;
					typeNum++;



					imgcopy[y*LineSize+x-1]=1;
				}
				//zuoshang
				if(x-1>=0&&y-1>=0&&p[(y-1)*LineSize+x-1]==DealPixel&&imgcopy[(y-1)*LineSize+x-1]==0)
				{
					typeImg[(y-1)*w+x-1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x-1;
					keytype[n].y=y-1;
					typeNum++;



					imgcopy[(y-1)*LineSize+x-1]=1;
				}
				//youshang
				if(y-1>=0&&x+1<w&&p[(y-1)*LineSize+x+1]==DealPixel&&imgcopy[(y-1)*LineSize+x+1]==0)
				{
					typeImg[(y-1)*w+x+1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x+1;
					keytype[n].y=y-1;
					typeNum++;



					imgcopy[(y-1)*LineSize+x+1]=1;
				}
				//you
				if(x+1<w&&p[y*LineSize+x+1]==DealPixel&&imgcopy[y*LineSize+x+1]==0)
				{
					typeImg[y*w+x+1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x+1;
					keytype[n].y=y;
					typeNum++;



					imgcopy[y*LineSize+x+1]=1;
				}
				//shang
				if(y-1>=0&&p[(y-1)*LineSize+x]==DealPixel&&imgcopy[(y-1)*LineSize+x]==0)
				{
					typeImg[(y-1)*w+x]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x;
					keytype[n].y=y-1;
					typeNum++;



					imgcopy[(y-1)*LineSize+x]=1;
				}
				//zuoxia
				if(x-1>=0 && y+1<h&&p[(y+1)*LineSize+x-1]==DealPixel&&imgcopy[(y+1)*LineSize+x-1]==0)
				{
					typeImg[(y+1)*w+x-1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x-1;
					keytype[n].y=y+1;
					typeNum++;



					imgcopy[(y+1)*LineSize+x-1]=1;
				}
				//xia
				if(y+1<h&&p[(y+1)*LineSize+x]==DealPixel&&imgcopy[(y+1)*LineSize+x]==0)
				{
					typeImg[(y+1)*w+x]=N+1;
					n++;


					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x;
					keytype[n].y=y+1;
					typeNum++;


					imgcopy[(y+1)*LineSize+x]=1;
				}
				//youxia
				if(y+1<h&&x+1<w&&p[(y+1)*LineSize+x+1]==DealPixel&&imgcopy[(y+1)*LineSize+x+1]==0)
				{
					typeImg[(y+1)*w+x+1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x+1;
					keytype[n].y=y+1;
					typeNum++;



					imgcopy[(y+1)*LineSize+x+1]=1;
				}
			}
			//过滤后保存在数组当中的裂纹为预选的裂纹,后续再进行合并
			if(typeNum>0)
			{
				if (N>=REGION_NUM)
				{
					free(imgcopy);
					return REGION_NUM;
				}

				N++;
				region[N].PointNum=typeNum;
				region[N].left=nLEFT;
				region[N].right=nRight;
				region[N].top=nTOP;
				region[N].bottom=nBOTTOM;
				region[N].IsOK=true;
				region[N].FusedN=0;
				region[N].W=nRight-nLEFT+1;
				region[N].H=nBOTTOM-nTOP+1;

				region[N].LeftTop=nTOP;
				region[N].RightTop=nTOP;
			}
		}
	}

	free(imgcopy);

	return ++N;
}

//双边滤波
//void BiFilter(IplImage * src, IplImage * des,int halfL,float delta1,float delta2)
//{
//	int i,j,m,n;
//	int height=src->height;
//	int width=src->width;
//	BYTE * Data=(BYTE *)src->imageData;
//
//	//高斯模糊的权重
//	int L=2*halfL+1;
//	float * gaussW=new float[L*L];
//	for (i=-halfL;i<=halfL;i++)
//	{
//		for (j=-halfL;j<=halfL;j++)
//		{
//			gaussW[(i+halfL)*L+(j+halfL)]=1/sqrt(2*3.1415926*delta1)*pow(2.718281828,-1.0*(i*i+j*j)/(2*delta1));
//		}
//	}
//
//	float * gaussPix=new float[256];
//	for (j=0;j<256;j++)
//	{
//		gaussPix[j]=1/sqrt(2*3.1415926*delta2)*pow(2.718281828,-1.0*j*j/(2*delta2));
//	}
//
//	//开始模糊
//	float wt,pixV,weightsum,temp;
//	for (i=0;i<height;i++)
//	{
//		for (j=0;j<width;j++)
//		{
//			for (int mm=0;mm<src->nChannels;mm++)
//			{
//				weightsum=0;
//				pixV=0;
//				for (m=-halfL;m<=halfL;m++)
//				{
//					for (n=-halfL;n<=halfL;n++)
//					{
//						if (i+m<0 || i+m>height-1 || j+n<0 ||j+n>width-1)
//						{
//
//						}
//						else
//						{
//							int pos=abs(Data[(i+m)*src->widthStep+(j+n)*src->nChannels+mm]-Data[i*src->widthStep+j*src->nChannels+mm]);
//							float W= (gaussW[(m+halfL)*L+(n+halfL)]*gaussPix[pos]);
//							weightsum+=W;
//							pixV=pixV+W*Data[(i+m)*src->widthStep+(j+n)*src->nChannels+mm];
//						}
//					}
//				}
//
//				des->imageData[i*des->widthStep+j*des->nChannels+mm]=(pixV+1e-6)/(weightsum+1e-6);
//			}
//		}
//	}
//
//	delete []gaussPix;
//	delete []gaussW;
//}
