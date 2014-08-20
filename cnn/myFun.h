#pragma once
int GetFileNameFromDir(char* _dir, char** filename, int filenum, char* suffix =NULL);

//初始化神经网络模型
void IniNN(NeuralNetwork& NN);

//提取样本集的单个样本
void ReadData(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, char* TestImages, char* TestLabels);

void ReadData(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, ifstream& fileTestImages, ifstream& fileTestLabels);


//该函数用来读取中文识别时的标签，因为中文类别超过2**8=256类别，需要用一个int来表明label
void ReadData_OCRZH(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, char* TestImages, char* TestLabels);

//载入分类器
void LoadNN(NeuralNetwork& NN, TCHAR* NN_name);

//进行识别，
void Recognition(NeuralNetwork& NN, unsigned char* pInput, double* pOutput, int oCount);

//扭曲输入
void DistortionMap(double *inputVector, double severityFactor = 1.0);


//根据样本文件夹制作开源Demo需要的样本文件
//调用此函数时注意修改文件对应标签在文件名中的位置（包括目录部分）
void MakeSamplesFile(char* SamplesDir, int SamplesNum, char* SamplesFileName);

//归一化
void normImg1(IplImage * src,IplImage * Des);
void normImg(IplImage * src,IplImage * Des);

void GetCharFromBackground(IplImage* kSrc, IplImage** pDst);