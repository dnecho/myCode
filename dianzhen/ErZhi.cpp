#include "stdafx.h"
#include "ErZhi.h"

#define S 23
#define T (0.15f)


void adaptiveThreshold(IplImage * GrayImg, IplImage * BiImg)
{
	if (GrayImg->nChannels>1)
	{
		cvCvtColor(GrayImg,GrayImg,CV_BGR2GRAY);
	}

	BYTE * input=(BYTE *)GrayImg->imageData;
	BYTE * bin=(BYTE *)BiImg->imageData;
	int IMAGE_WIDTH=GrayImg->width;
	int IMAGE_HEIGHT=GrayImg->height;
	int IMAGE_WIDESTEP=GrayImg->widthStep;

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

void SML_OTSU(IplImage * Src,IplImage * Des)
{
	if (Src==NULL || Des==NULL)
	{
		return;
	}

	int i,j;

	IplImage * GrayImg=NULL;
	if (Src->nChannels>1)
	{
		GrayImg=cvCreateImage(cvGetSize(Src),8,1);
		cvCvtColor(Src,GrayImg,CV_BGR2GRAY);
	}
	else
	{
		GrayImg=cvCreateImage(cvGetSize(Src),8,1);
		cvCopy(Src,GrayImg);
	}

	int height=GrayImg->height;
	int width=GrayImg->width;
	int pixelSum=height*width;
	int widestep=GrayImg->widthStep;
	BYTE * Data=(BYTE *)GrayImg->imageData;

	int pixelCount[256];  
	float pixelPro[256];  

	memset(pixelCount,0,sizeof(int)*256);
	memset(pixelPro,0,sizeof(float)*256);

	//统计灰度级中每个像素在整幅图像中的个数  
	for (i=0;i<height;i++)
	{
		for (j=0;j<width;j++)
		{
			pixelCount[(int)Data[i*widestep+j]]++;
		}
	}

	for (i=0;i<256;i++)
	{
		pixelPro[i]=1.0f*pixelCount[i]/pixelSum;
	}

	float w0, w1, u0tmp, u1tmp, u0, u1, u, deltaTmp, deltaMax = 0;  
	int threshold;
	for (i=0;i<256;i++)
	{
		w0 = w1 = u0tmp = u1tmp = u0 = u1 = u = deltaTmp = 0;  
		for (j=0;j<256;j++)
		{
			if (j<=i)
			{
				w0+=pixelPro[j];
				u0tmp+=j*pixelPro[j];
			}
			else
			{
				w1+=pixelPro[j];
				u1tmp+=j*pixelPro[j];
			}
		}

		if (w0==0)
			u0=0;
		else
			u0 = u0tmp / w0; 
		if (w1==0)
			u1=0;
		else
			u1 = u1tmp / w1;  
		u = u0tmp + u1tmp;  
		deltaTmp = w0 * pow((u0 - u), 2) + w1 * pow((u1 - u), 2);  
		if(deltaTmp > deltaMax)  
		{  
			deltaMax = deltaTmp;  
			threshold = i;  
		}  
	}

	//经验值调整
	threshold-=22;

	BYTE * sData=(BYTE *)Des->imageData;
	for (i=0;i<height;i++)
	{
		for (j=0;j<width;j++)
		{
			if (Data[i*widestep+j]>=threshold)
			{
				sData[i*widestep+j]=255;
			}
			else
			{
				sData[i*widestep+j]=0;
			}
		}
	}
	cvReleaseImage(&GrayImg);
}


void SML_OTSU_Bersen(IplImage * Src,IplImage * Des)
{
	if (Src==NULL || Des==NULL)
	{
		return;
	}

	int i,j,m,n;

	IplImage * GrayImg=NULL;
	if (Src->nChannels>1)
	{
		GrayImg=cvCreateImage(cvGetSize(Src),8,1);
		cvCvtColor(Src,GrayImg,CV_BGR2GRAY);
	}
	else
	{
		GrayImg=cvCreateImage(cvGetSize(Src),8,1);
		cvCopy(Src,GrayImg);
	}

	int height=GrayImg->height;
	int width=GrayImg->width;
	int pixelSum=height*width;
	int widestep=GrayImg->widthStep;
	BYTE * Data=(BYTE *)GrayImg->imageData;

	int pixelCount[256];  
	float pixelPro[256];  

	memset(pixelCount,0,sizeof(int)*256);
	memset(pixelPro,0,sizeof(float)*256);

	//统计灰度级中每个像素在整幅图像中的个数  
	for (i=0;i<height;i++)
	{
		for (j=0;j<width;j++)
		{
			pixelCount[(int)Data[i*widestep+j]]++;
		}
	}

	for (i=0;i<256;i++)
	{
		pixelPro[i]=1.0f*pixelCount[i]/pixelSum;
	}

	float w0, w1, u0tmp, u1tmp, u0, u1, u, deltaTmp, deltaMax = 0;  
	int threshold;
	for (i=0;i<256;i++)
	{
		w0 = w1 = u0tmp = u1tmp = u0 = u1 = u = deltaTmp = 0;  
		for (j=0;j<256;j++)
		{
			if (j<=i)
			{
				w0+=pixelPro[j];
				u0tmp+=j*pixelPro[j];
			}
			else
			{
				w1+=pixelPro[j];
				u1tmp+=j*pixelPro[j];
			}
		}

		if (w0==0)
			u0=0;
		else
			u0 = u0tmp / w0; 
		if (w1==0)
			u1=0;
		else
			u1 = u1tmp / w1;  
		u = u0tmp + u1tmp;  
		deltaTmp = w0 * pow((u0 - u), 2) + w1 * pow((u1 - u), 2);  
		if(deltaTmp > deltaMax)  
		{  
			deltaMax = deltaTmp;  
			threshold = i;  
		}  
	}

	int WindowRadio=21;
	int r=WindowRadio/2;
	int TH=20;
	memset(Des->imageData,255,Des->height*Des->widthStep);
	BYTE * sData=(BYTE *)Des->imageData;
	for (i=r;i<height-r;i++)
	{
		for (j=r;j<width-r;j++)
		{
			//将距离较远的点，根据全局阈值进行分割
			if (abs(Data[i*widestep+j]-threshold)>TH)
			{
				if (Data[i*widestep+j]>threshold)
				{
					sData[i*widestep+j]=255;
				}
				else
					sData[i*widestep+j]=0;
			}
			//距离相近的点，根据周围点的信息进行分割
			else
			{
				BYTE B=Data[i*widestep+j*1];
				int minP=255;
				int maxP=0;
				for (m=i-r;m<=i+r;m++)
				{
					for (n=j-r;n<=j+r;n++)
					{
						if (Data[m*widestep+n]<minP)
						{
							minP=Data[m*widestep+n];
						}

						if (Data[m*widestep+n]>maxP)
						{
							maxP=Data[m*widestep+n];
						}
					}
				}

				if (maxP-B>B-minP)
				{
					sData[i*Des->widthStep+j]=0;
				}
				else
				{
					sData[i*Des->widthStep+j]=255;
				}
			}
		}
	}

	cvReleaseImage(&GrayImg);
}