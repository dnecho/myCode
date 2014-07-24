#ifndef ROTATE_H
#define ROTATE_H

#define PI 3.1415926535
#define RADIAN(angle) ((angle)*PI/180.0) //�Ƕȵ�����ת���ĺ�
typedef unsigned char BYTE;

IplImage * RotateGrayImg(IplImage * src, double angle);

IplImage* RotateImage_CV(IplImage* img,float degree);
#endif