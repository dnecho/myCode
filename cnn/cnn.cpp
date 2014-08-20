#include "GlobalVal.h"


int _tmain(int argc, _TCHAR* argv[])
{

	NeuralNetwork myTestCNN;
	myTestCNN.LoadNN_C("lic_C.nnt");

	//��ȡ����
	int idx = 1;
	int label = -1;
	char* TestImages = "DigPrintImages4443.idx3-ubyte";
	char* TestLabels = "DigPrintLabels4443.idx1-ubyte";

	//��������
	ifstream fileTestImages;
	fileTestImages.open(TestImages, ios::in|ios::binary);

	//��������ʵ�ʱ�ǩ
	ifstream fileTestLabels;
	fileTestLabels.open(TestLabels, ios::in|ios::binary);

	int SampleNum = 4443;
	int WrongNum = 0;
	unsigned char pArray[ g_cImageSize * g_cImageSize ] = {0};
	for(int i = 0; i < SampleNum; i++)
	{
		cout<<i<<":";
		ReadData(i, pArray, &label, true, fileTestImages, fileTestLabels);


		//ʶ������
		//DWORD tick = ::GetTickCount();
		const int outNode = 10;//����ڵ����
		double outputVector[outNode] = {0.0};
		//Recognition(myTestCNN, pArray, outputVector, outNode);
		myTestCNN.Recognition(pArray, outputVector, outNode);
		//cout<<"time:"<<GetTickCount()-tick<<endl;

		int maxpos;
		double maxV=-10000;
		for(int i = 0; i < outNode; i++)
		{
			//cout<<outputVector[i]<<"\t";

			if (outputVector[i]>maxV)
			{
				maxV=outputVector[i];
				maxpos=i;
			}
		}
		//cout<<endl;

		//int pow_2[12]={2048,1024,512,256,128,64,32,16,8,4,2,1};
		//maxpos = 0; 
		//for(int ii=0; ii<12; ++ii)
		//{
		//	if(outputVector[ii] > 0.5)
		//	{
		//		maxpos += pow_2[ii];
		//	}
		//}

		cout<<"Label:"<<label<<"\t";
		cout<<"Result:"<<maxpos<<"\t";
		if(maxpos != label)
		{
			cout<<"����\t"<<endl;
			WrongNum++;
		}
		else
			cout<<"\t��ȷ"<<endl;
	}
	cout<<"����"<<WrongNum<<"��,ʶ������"<<1 - WrongNum*0.1/SampleNum<<endl;


	system("pause");
	return 0;
}
