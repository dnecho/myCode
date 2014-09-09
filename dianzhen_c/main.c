#include "dianzhen.h"
#include "region.h"
#include <stdio.h>
#include <time.h>

int main(){
	
	FILE *fp = fopen(".//pictxt//qvga1.txt","r");
	int w = 0;
	int h = 0;
	int i,j,tmp;
	int IsSuccess;
	unsigned char* pImgData;
	Axis center={0};
	int pixelsize=0;
	long t_start,t_end;
	
	fscanf(fp, "%d", &w);
	fscanf(fp, "%d", &h);
	
	pImgData = (unsigned char*)malloc(sizeof(unsigned char)*w*h);
	
	
	for(i = 0; i < h; i++)
	{
		for(j = 0; j < w; j++)
		{
			fscanf(fp, "%d", &tmp);
			pImgData[j+w*i]=(unsigned char)tmp;
		}
	}
	
	t_start = clock();
	IsSuccess = process(pImgData, w, h, &center, &pixelsize);
	if(IsSuccess)
	{
		printf("final\tx:%d£¬y:%d£¬PixelSize:%d\n", center.x, center.y, pixelsize);
	}else
	{
		printf("´¦ÀíÊ§°Ü\n");
	}
	t_end = clock();
	//printf("total:%d", t_end-t_start);
	//getchar();
	system("pause");

}