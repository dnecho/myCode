#pragma once
int GetFileNameFromDir(char* _dir, char** filename, int filenum, char* suffix =NULL);

//��ʼ��������ģ��
void IniNN(NeuralNetwork& NN);

//��ȡ�������ĵ�������
void ReadData(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, char* TestImages, char* TestLabels);

void ReadData(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, ifstream& fileTestImages, ifstream& fileTestLabels);


//�ú���������ȡ����ʶ��ʱ�ı�ǩ����Ϊ������𳬹�2**8=256�����Ҫ��һ��int������label
void ReadData_OCRZH(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, char* TestImages, char* TestLabels);

//���������
void LoadNN(NeuralNetwork& NN, TCHAR* NN_name);

//����ʶ��
void Recognition(NeuralNetwork& NN, unsigned char* pInput, double* pOutput, int oCount);

//Ť������
void DistortionMap(double *inputVector, double severityFactor = 1.0);


//���������ļ���������ԴDemo��Ҫ�������ļ�
//���ô˺���ʱע���޸��ļ���Ӧ��ǩ���ļ����е�λ�ã�����Ŀ¼���֣�
void MakeSamplesFile(char* SamplesDir, int SamplesNum, char* SamplesFileName);

//��һ��
void normImg1(IplImage * src,IplImage * Des);
void normImg(IplImage * src,IplImage * Des);

void GetCharFromBackground(IplImage* kSrc, IplImage** pDst);