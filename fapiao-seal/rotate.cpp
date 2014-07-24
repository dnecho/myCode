#include "stdafx.h"
IplImage * RotateGrayImg(IplImage * src, double angle)
{
//	int m_Angle=angle;
	int width=src->width;
	int height=src->height;

	BYTE * lpTempPtr=NULL;
	BYTE * p_Copy=(BYTE *)src->imageData;

	double       SrcX1,SrcY1,SrcX2,SrcY2;
	double       SrcX3,SrcY3,SrcX4,SrcY4;
	double       DstX1,DstY1,DstX2,DstY2;
	double       DstX3,DstY3,DstX4,DstY4;
	double       num1,num2;

	double RotateAngle=(double)RADIAN(angle);
	//double RotateAngle=angle;

	double cosa=(double)cos((double)RotateAngle);
	double sina=(double)sin((double)RotateAngle);
	//ԭͼ�Ŀ�Ⱥ͸߶�
	//   BYTE * lpSrc=p_Copy;
	int Wold=width;
	int Hold=height;
	//ԭͼ���ĸ��ǵ�����
	SrcX1=(double)(-0.5*Wold);
	SrcY1=(double)(0.5*Hold);
	SrcX2=(double)(0.5*Wold);
	SrcY2=(double)(0.5*Hold);
	SrcX3=(double)(-0.5*Wold);
	SrcY3=(double)(-0.5*Hold);
	SrcX4=(double)(0.5*Wold);
	SrcY4=(double)(-0.5*Hold);
	//��ͼ�ĸ��ǵ�����
	DstX1=cosa*SrcX1+sina*SrcY1;
	DstY1=-sina*SrcX1+cosa*SrcY1;
	DstX2=cosa*SrcX2+sina*SrcY2;
	DstY2=-sina*SrcX2+cosa*SrcY2;
	DstX3=cosa*SrcX3+sina*SrcY3;
	DstY3=-sina*SrcX3+cosa*SrcY3;
	DstX4=cosa*SrcX4+sina*SrcY4;
	DstY4=-sina*SrcX4+cosa*SrcY4;
	//������ͼ�Ŀ�ȣ��߶�
	int Wnew = (int)(max(fabs(DstX4-DstX1), fabs(DstX3-DstX2))+0.5);
	int Hnew = (int)(max(fabs(DstY4-DstY1), fabs(DstY3-DstY2))+0.5);
	//�������(2.9)�е��������������������Ժ�ÿ�ζ�������
	num1=(double)( -0.5*Wnew*cosa-0.5*Hnew*sina+0.5*Wold);
	num2=(double)(0.5*Wnew*sina-0.5*Hnew*cosa+0.5*Hold);

	//���µĻ������е�ÿ���ֽڶ����255�������Ժ�δ��������ؾ��ǰ�ɫ
	int DstBufSize=Wnew*Hnew;

	/*if(lpTempPtr==NULL)
	{
		lpTempPtr=new BYTE[Wnew*Hnew];
		memset(lpTempPtr,(BYTE)255,Wnew*Hnew);
	}
	else
	{
		delete []lpTempPtr;
		lpTempPtr=NULL;

		lpTempPtr=new BYTE[Wnew*Hnew];
		memset(lpTempPtr,(BYTE)255,Wnew*Hnew);
	}*/

	IplImage *Des=cvCreateImage(cvSize(Wnew,Hnew),src->depth,src->nChannels);
	lpTempPtr=(BYTE *)Des->imageData;
	memset(lpTempPtr,255,Des->widthStep*Des->height);

	int        x0,y0,x1,y1;
	int nChannels = src->nChannels;
	for(y1=0;y1<Hnew;y1++)
	{
		for(x1=0;x1<Wnew;x1++)
		{
			//x0,y0Ϊ��Ӧ��ԭͼ�ϵ�����
			//for(int iChannel = 0; iChannel < nChannels; iChannel++)
			//{
			//	x0= (int)((x1*nChannels+iChannel)*cosa+y1*sina+num1);
			//	y0= (int)(-1.0f*(x1*nChannels+iChannel)*sina+y1*cosa+num2);
			//	if( (x0>=0) && (x0<src->widthStep) && (y0>=0) && (y0<Hold))   //��ԭͼ��Χ��
			//	{
			//		lpTempPtr[y1*Des->widthStep+x1*nChannels+iChannel]=p_Copy[y0*src->widthStep+x0];
			//	}
			//}

			x0= (int)(x1*cosa+y1*sina+num1);
			y0= (int)(-1.0f*x1*sina+y1*cosa+num2);
			if( (x0>=0) && (x0<Wold) && (y0>=0) && (y0<Hold))   //��ԭͼ��Χ��
			{
				lpTempPtr[y1*Des->widthStep+x1]=p_Copy[y0*src->widthStep+x0];
			}
		}
	}

//	MyNamedWindow("rr",1);
//	MyShowImage("rr",Des);

	return Des;
}

IplImage* RotateImage_CV(IplImage* img,float degree)
{  
	double angle = degree  * CV_PI / 180.; // ����    
	double a = sin(angle), b = cos(angle);   
	int width = img->width;    
	int height = img->height;    
	int width_rotate= int(height * fabs(a) + width * fabs(b));    
	int height_rotate=int(width * fabs(a) + height * fabs(b));    
	//��ת����map  
	// [ m0  m1  m2 ] ===>  [ A11  A12   b1 ]  
	// [ m3  m4  m5 ] ===>  [ A21  A22   b2 ]  
	float map[6]={0};  
	CvMat map_matrix = cvMat(2, 3, CV_32F, map);    
	// ��ת����  
	CvPoint2D32f center = cvPoint2D32f(width / 2, height / 2);    
	cv2DRotationMatrix(center, degree, 1.0, &map_matrix);    
	map[2] += (width_rotate - width) / 2;    
	map[5] += (height_rotate - height) / 2;    
	IplImage* img_rotate = cvCreateImage(cvSize(width_rotate, height_rotate), 8, img->nChannels);   
	//��ͼ��������任  
	//CV_WARP_FILL_OUTLIERS - ����������ͼ������ء�  
	//�������������������ͼ��ı߽��⣬��ô���ǵ�ֵ�趨Ϊ fillval.  
	//CV_WARP_INVERSE_MAP - ָ�� map_matrix �����ͼ������ͼ��ķ��任��  
	cvWarpAffine( img,img_rotate, &map_matrix, CV_INTER_LINEAR | CV_WARP_FILL_OUTLIERS, cvScalarAll(0));    
	return img_rotate;  
}  