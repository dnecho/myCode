#include "region.h"
#include <memory.h>
#include <malloc.h>
Region region[1000];
#define MaxPointNum 20000
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

	//���ݾ�������ж��Ƿ���Ҫ�������ֽڶ���
	LineSize=(w+3)/4*4;
	DealPixel=0;

	imgcopy=(int*)malloc(sizeof(int)*LineSize*h);
	memset(imgcopy,0,sizeof(int)*LineSize*h);

	//const int MaxPointNum=20000;
	//CPoint keytype[MaxPointNum];

	//����С��
	N=-1;
	//int nLEFT,nRight,nTOP,nBOTTOM;
	for (j=0;j<h;j++)
	{
		for (i=0;i<w;i++)
		{
			//�����±�
			n=-1;
			//��ǰ�׵����Ŀ
			typeNum=0;
			//�����ǰ���ǰ׵㣬���һ�û�д������������Ϊһ�����ӵ㣬��������ɢ
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
				//����ĩβ��ֵ�൱��POP����
				x=keytype[n].x;
				y=keytype[n].y;

				nLEFT=nLEFT<x?nLEFT:x;
				nRight=nRight>x?nRight:x;
				nTOP=nTOP<y?nTOP:y;
				nBOTTOM=nBOTTOM>y?nBOTTOM:y;

				//�±��һ
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
			//���˺󱣴������鵱�е�����ΪԤѡ������,�����ٽ��кϲ�
			if(typeNum>0)
			{
				if (N>=1000)
				{
					free(imgcopy);
					return 1000;
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
