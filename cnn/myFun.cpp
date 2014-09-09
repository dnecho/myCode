#include "GlobalVal.h"

//dir:文件目录最后不要求为\\，如"E:\\face"
//filename:文件名字，包含目录，可直接用于打开文件
//filenum:希望从目录中得到的文件数量
//suffix:希望得到文件的后缀，如jpg，bmp（注意不用加.），默认为空，为空时不指定文件后缀名
//函数返回值表示实际得到的文件数量
int GetFileNameFromDir(char* _dir, char** filename, int filenum, char* suffix)
{
	int num = 0;
	int len=strlen(_dir);
	char dir[255];
	if(dir[len - 1] != '\\')
		sprintf_s(dir, "%s\\",_dir);
	else
		strcpy_s(dir, _dir);
	_finddata_t file;
	long longf;
	char tep_char[255];

	//得到需要查找的路径及需要查找文件的类型（是否指定后缀名）
	if(suffix == NULL)//判断是否需要制定文件后缀
	{
		sprintf_s(tep_char, "%s*.*", dir);
	}
	else
	{		
		sprintf_s(tep_char, "%s*.%s", dir, suffix);
	}

	//_findfirst返回的是long型; long __cdecl _findfirst(const char *, struct _finddata_t *)
	if((longf = _findfirst(tep_char, &file))==-1l)
	{
		printf("文件没有找到!/n");
	}
	else
	{
		if (file.name[0] != '.')//当后缀是"null"时findfirst找到的是'.'；如果后缀有值是找到的是第一个对应后缀的文件，必须保存起来
		{
			sprintf_s(tep_char, "%s%s", dir, file.name);//将路径与文件名结合
			strcpy_s(filename[num], 255, tep_char);
			num++;
		}

		//int __cdecl _findnext(long, struct _finddata_t *);如果找到下个文件的名字成功的话就返回0,否则返回-1
		while( _findnext( longf, &file ) == 0 && num <filenum)
		{
			if (file.name[0] == '.' && file.name[1] == '.')
			{
				continue;
			}
			sprintf_s(tep_char, "%s%s", dir, file.name);//将路径与文件名结合
			strcpy_s(filename[num], 255, tep_char);
			num++;
		}
	}

	_findclose(longf);
	return num;
}


//void IniNN(NeuralNetwork& NN)
//{
//
//	// initialize and build the neural net
//
//	NN.Initialize();
//
//	NNLayer* pLayer;
//
//	int ii, jj, kk;
//	int icNeurons = 0;
//	int icWeights = 0;
//	double initWeight;
//	string label;
//
//	// layer zero, the input layer.
//	// Create neurons: exactly the same number of neurons as the input
//	// vector of 29x29=841 pixels, and no weights/connections
//
//	pLayer = new NNLayer( _T("Layer00") );
//	NN.m_Layers.push_back( pLayer );
//
//	for ( ii=0; ii<841; ++ii )
//	{
//		label.Format( _T("Layer00_Neuron%04d_Num%06d"), ii, icNeurons );
//		pLayer->m_Neurons.push_back( new NNNeuron( (LPCTSTR)label ) );
//		icNeurons++;
//	}
//
//#define UNIFORM_PLUS_MINUS_ONE ( (double)(2.0 * rand())/RAND_MAX - 1.0 )
//
//	// layer one:
//	// This layer is a convolutional layer that has 6 feature maps.  Each feature 
//	// map is 13x13, and each unit in the feature maps is a 5x5 convolutional kernel
//	// of the input layer.
//	// So, there are 13x13x6 = 1014 neurons, (5x5+1)x6 = 156 weights
//
//	pLayer = new NNLayer( _T("Layer01"), pLayer );
//	NN.m_Layers.push_back( pLayer );
//
//	for ( ii=0; ii<1014; ++ii )
//	{
//		label.Format( _T("Layer01_Neuron%04d_Num%06d"), ii, icNeurons );
//		pLayer->m_Neurons.push_back( new NNNeuron( (LPCTSTR)label ) );
//		icNeurons++;
//	}
//
//	for ( ii=0; ii<156; ++ii )
//	{
//		label.Format( _T("Layer01_Weight%04d_Num%06d"), ii, icWeights );
//		initWeight = 0.05 * UNIFORM_PLUS_MINUS_ONE;
//		pLayer->m_Weights.push_back( new NNWeight( (LPCTSTR)label, initWeight ) );
//	}
//
//	// interconnections with previous layer: this is difficult
//	// The previous layer is a top-down bitmap image that has been padded to size 29x29
//	// Each neuron in this layer is connected to a 5x5 kernel in its feature map, which 
//	// is also a top-down bitmap of size 13x13.  We move the kernel by TWO pixels, i.e., we
//	// skip every other pixel in the input image
//
//	int kernelTemplate[25] = {
//		0,  1,  2,  3,  4,
//		29, 30, 31, 32, 33,
//		58, 59, 60, 61, 62,
//		87, 88, 89, 90, 91,
//		116,117,118,119,120 };
//
//		int iNumWeight;
//
//		int fm;
//
//		for ( fm=0; fm<6; ++fm)
//		{
//			for ( ii=0; ii<13; ++ii )
//			{
//				for ( jj=0; jj<13; ++jj )
//				{
//					iNumWeight = fm * 26;  // 26 is the number of weights per feature map
//					NNNeuron& n = *( pLayer->m_Neurons[ jj + ii*13 + fm*169 ] );
//
//					n.AddConnection( ULONG_MAX, iNumWeight++ );  // bias weight
//
//					for ( kk=0; kk<25; ++kk )
//					{
//						// note: max val of index == 840, corresponding to 841 neurons in prev layer
//						n.AddConnection( 2*jj + 58*ii + kernelTemplate[kk], iNumWeight++ );
//					}
//				}
//			}
//		}
//
//
//		// layer two:
//		// This layer is a convolutional layer that has 50 feature maps.  Each feature 
//		// map is 5x5, and each unit in the feature maps is a 5x5 convolutional kernel
//		// of corresponding areas of all 6 of the previous layers, each of which is a 13x13 feature map
//		// So, there are 5x5x50 = 1250 neurons, (5x5+1)x6x50 = 7800 weights
//
//		pLayer = new NNLayer( _T("Layer02"), pLayer );
//		NN.m_Layers.push_back( pLayer );
//
//		for ( ii=0; ii<1250; ++ii )
//		{
//			label.Format( _T("Layer02_Neuron%04d_Num%06d"), ii, icNeurons );
//			pLayer->m_Neurons.push_back( new NNNeuron( (LPCTSTR)label ) );
//			icNeurons++;
//		}
//
//		for ( ii=0; ii<7800; ++ii )
//		{
//			label.Format( _T("Layer02_Weight%04d_Num%06d"), ii, icWeights );
//			initWeight = 0.05 * UNIFORM_PLUS_MINUS_ONE;
//			pLayer->m_Weights.push_back( new NNWeight( (LPCTSTR)label, initWeight ) );
//		}
//
//		// Interconnections with previous layer: this is difficult
//		// Each feature map in the previous layer is a top-down bitmap image whose size
//		// is 13x13, and there are 6 such feature maps.  Each neuron in one 5x5 feature map of this 
//		// layer is connected to a 5x5 kernel positioned correspondingly in all 6 parent
//		// feature maps, and there are individual weights for the six different 5x5 kernels.  As
//		// before, we move the kernel by TWO pixels, i.e., we
//		// skip every other pixel in the input image.  The result is 50 different 5x5 top-down bitmap
//		// feature maps
//
//		int kernelTemplate2[25] = {
//			0,  1,  2,  3,  4,
//			13, 14, 15, 16, 17, 
//			26, 27, 28, 29, 30,
//			39, 40, 41, 42, 43, 
//			52, 53, 54, 55, 56   };
//
//
//			for ( fm=0; fm<50; ++fm)
//			{
//				for ( ii=0; ii<5; ++ii )
//				{
//					for ( jj=0; jj<5; ++jj )
//					{
//						iNumWeight = fm * 26;  // 26 is the number of weights per feature map
//						NNNeuron& n = *( pLayer->m_Neurons[ jj + ii*5 + fm*25 ] );
//
//						n.AddConnection( ULONG_MAX, iNumWeight++ );  // bias weight
//
//						for ( kk=0; kk<25; ++kk )
//						{
//							// note: max val of index == 1013, corresponding to 1014 neurons in prev layer
//							n.AddConnection(       2*jj + 26*ii + kernelTemplate2[kk], iNumWeight++ );
//							n.AddConnection( 169 + 2*jj + 26*ii + kernelTemplate2[kk], iNumWeight++ );
//							n.AddConnection( 338 + 2*jj + 26*ii + kernelTemplate2[kk], iNumWeight++ );
//							n.AddConnection( 507 + 2*jj + 26*ii + kernelTemplate2[kk], iNumWeight++ );
//							n.AddConnection( 676 + 2*jj + 26*ii + kernelTemplate2[kk], iNumWeight++ );
//							n.AddConnection( 845 + 2*jj + 26*ii + kernelTemplate2[kk], iNumWeight++ );
//						}
//					}
//				}
//			}
//
//
//			// layer three:
//			// This layer is a fully-connected layer with 100 units.  Since it is fully-connected,
//			// each of the 100 neurons in the layer is connected to all 1250 neurons in
//			// the previous layer.
//			// So, there are 100 neurons and 100*(1250+1)=125100 weights
//
//			pLayer = new NNLayer( _T("Layer03"), pLayer );
//			NN.m_Layers.push_back( pLayer );
//
//			for ( ii=0; ii<100; ++ii )
//			{
//				label.Format( _T("Layer03_Neuron%04d_Num%06d"), ii, icNeurons );
//				pLayer->m_Neurons.push_back( new NNNeuron( (LPCTSTR)label ) );
//				icNeurons++;
//			}
//
//			for ( ii=0; ii<125100; ++ii )
//			{
//				label.Format( _T("Layer03_Weight%04d_Num%06d"), ii, icWeights );
//				initWeight = 0.05 * UNIFORM_PLUS_MINUS_ONE;
//				pLayer->m_Weights.push_back( new NNWeight( (LPCTSTR)label, initWeight ) );
//			}
//
//			// Interconnections with previous layer: fully-connected
//
//			iNumWeight = 0;  // weights are not shared in this layer
//
//			for ( fm=0; fm<100; ++fm )
//			{
//				NNNeuron& n = *( pLayer->m_Neurons[ fm ] );
//				n.AddConnection( ULONG_MAX, iNumWeight++ );  // bias weight
//
//				for ( ii=0; ii<1250; ++ii )
//				{
//					n.AddConnection( ii, iNumWeight++ );
//				}
//			}
//
//
//
//			// layer four, the final (output) layer:
//			// This layer is a fully-connected layer with 10 units.  Since it is fully-connected,
//			// each of the 10 neurons in the layer is connected to all 100 neurons in
//			// the previous layer.
//			// So, there are 10 neurons and 10*(100+1)=1010 weights
//
//			pLayer = new NNLayer( _T("Layer04"), pLayer );
//			NN.m_Layers.push_back( pLayer );
//
//			for ( ii=0; ii<10; ++ii )
//			{
//				label.Format( _T("Layer04_Neuron%04d_Num%06d"), ii, icNeurons );
//				pLayer->m_Neurons.push_back( new NNNeuron( (LPCTSTR)label ) );
//				icNeurons++;
//			}
//
//			for ( ii=0; ii<1010; ++ii )
//			{
//				label.Format( _T("Layer04_Weight%04d_Num%06d"), ii, icWeights );
//				initWeight = 0.05 * UNIFORM_PLUS_MINUS_ONE;
//				pLayer->m_Weights.push_back( new NNWeight( (LPCTSTR)label, initWeight ) );
//			}
//
//			// Interconnections with previous layer: fully-connected
//
//			iNumWeight = 0;  // weights are not shared in this layer
//
//			for ( fm=0; fm<10; ++fm )
//			{
//				NNNeuron& n = *( pLayer->m_Neurons[ fm ] );
//				n.AddConnection( ULONG_MAX, iNumWeight++ );  // bias weight
//
//				for ( ii=0; ii<100; ++ii )
//				{
//					n.AddConnection( ii, iNumWeight++ );
//				}
//			}
//}
//

void ReadData(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, char* TestImages, char* TestLabels)
{
	int cCount = g_cImageSize*g_cImageSize;//总的数据
	int fPos;//偏移量

	//测试数据
	ifstream fileTestImages;
	fileTestImages.open(TestImages, ios::in|ios::binary);

	//测试数据实际标签
	ifstream fileTestLabels;
	fileTestLabels.open(TestLabels, ios::in|ios::binary);

	if ( pArray != NULL )
	{
		fPos = 16 + iNumImage*cCount;  // 16 compensates for file header info
		fileTestImages.seekg( fPos, ios::beg );
		fileTestImages.read( (char*)pArray, cCount );

		if ( bFlipGrayscale != false )
		{
			for ( int ii=0; ii<cCount; ++ii )
			{
				pArray[ ii ] = 255 - pArray[ ii ];
			}
		}
	}

	if ( pLabel != NULL )
	{
		fPos = 8 + iNumImage;
		char r;
		fileTestLabels.seekg( fPos, ios::beg );
		fileTestLabels.read( &r, 1 );  // single byte

		*pLabel = r;
	}	
}

void ReadData(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, ifstream& fileTestImages, ifstream& fileTestLabels)
{
	int cCount = g_cImageSize*g_cImageSize;//总的数据
	int fPos;//偏移量

	////测试数据
	//ifstream fileTestImages;
	//fileTestImages.open(TestImages, ios::in|ios::binary);

	////测试数据实际标签
	//ifstream fileTestLabels;
	//fileTestLabels.open(TestLabels, ios::in|ios::binary);

	if ( pArray != NULL )
	{
		fPos = 16 + iNumImage*cCount;  // 16 compensates for file header info
		fileTestImages.seekg( fPos, ios::beg );
		fileTestImages.read( (char*)pArray, cCount );

		if ( bFlipGrayscale != false )
		{
			for ( int ii=0; ii<cCount; ++ii )
			{
				pArray[ ii ] = 255 - pArray[ ii ];
			}
		}
	}

	if ( pLabel != NULL )
	{
		fPos = 8 + iNumImage;
		char r;
		fileTestLabels.seekg( fPos, ios::beg );
		fileTestLabels.read( &r, 1 );  // single byte

		*pLabel = r;
	}	
}

void ReadData_OCRZH(int iNumImage /* =0 */, unsigned char *pArray /* =NULL */, int *pLabel /* =NULL */,
	bool bFlipGrayscale /* =TRUE */, char* TestImages, char* TestLabels)
{
	int cCount = g_cImageSize*g_cImageSize;//总的数据
	int fPos;//偏移量

	//测试数据
	ifstream fileTestImages;
	fileTestImages.open(TestImages, ios::in|ios::binary);

	//测试数据实际标签
	ifstream fileTestLabels;
	fileTestLabels.open(TestLabels, ios::in|ios::binary);

	if ( pArray != NULL )
	{
		fPos = 16 + iNumImage*cCount;  // 16 compensates for file header info
		fileTestImages.seekg( fPos, ios::beg );
		fileTestImages.read( (char*)pArray, cCount );

		if ( bFlipGrayscale != false )
		{
			for ( int ii=0; ii<cCount; ++ii )
			{
				pArray[ ii ] = 255 - pArray[ ii ];
			}
		}
	}

	if ( pLabel != NULL )
	{
		fPos = 8 + iNumImage * sizeof(int);
		int r;
		fileTestLabels.seekg( fPos, ios::beg );
		fileTestLabels.read( (char*)&r, sizeof(int) );  // single byte

		*pLabel = r;
	}	
}
//
//void LoadNN(NeuralNetwork& NN, TCHAR* NN_name)
//{
//	CFile fileCNN;
//	if(fileCNN.Open(NN_name,CFile::modeRead)==0)
//	{
//		printf("cann't open classfier file(%s)", NN_name);
//		return;
//	}
//	CArchive ar(&fileCNN, CArchive::load);
//	NN.Initialize();//删除原有的神经网络模型
//	NN.Serialize(ar);
//}

void Recognition(NeuralNetwork& NN, unsigned char* pInput, double* pOutput, int oCount)
{
	int ii, jj;
	double inputVector[841];

	for ( ii=0; ii<841; ++ii )
	{
		inputVector[ ii ] = 1.0;  // one is white, -one is black
	}

	// top row of inputVector is left as zero, left-most column is left as zero 

	for ( ii=0; ii<g_cImageSize; ++ii )
	{
		for ( jj=0; jj<g_cImageSize; ++jj )
		{
			inputVector[ 1 + jj + 29*(ii+1) ] = (double)((int)(unsigned char)pInput[ jj + g_cImageSize*ii ])/128.0 - 1.0;  // one is white, -one is black
		}
	}
	
	NN.Calculate(inputVector, 841, pOutput, oCount, NULL);
}


//inline DWORD& At( DWORD* p, int row, int col )  // zero-based indices, starting at bottom-left
//{	
//	int location = row * m_cCols + col;
//	ASSERT( location>=0 && location<m_cPixels && row<m_cRows && row>=0 && col<m_cCols && col>=0 );
//	return p[ location ];
//}

inline double& At( double* p, int row, int col )  // zero-based indices, starting at bottom-left
{ 
	int location = row * m_cCols + col;
	ASSERT( location>=0 && location<m_cCount && row<m_cRows && row>=0 && col<m_cCols && col>=0 );
	return p[ location ];
}

//函数未完成，有很多成员变量需要处理，后期有需要再进行修改
//void DistortionMap(double severityFactor, double *inputVector)
//{
//	// generates distortion maps in each of the horizontal and vertical directions
//	// Three distortions are applied: a scaling, a rotation, and an elastic distortion
//	// Since these are all linear tranformations, we can simply add them together, after calculation
//	// one at a time
//
//	// The input parameter, severityFactor, let's us control the severity of the distortions relative
//	// to the default values.  For example, if we only want half as harsh a distortion, set
//	// severityFactor == 0.5
//
//	// First, elastic distortion, per Patrice Simard, "Best Practices For Convolutional Neural Networks..."
//	// at page 2.
//	// Three-step process: seed array with uniform randoms, filter with a gaussian kernel, normalize (scale)
//
//	int row, col;
//	double* uniformH = new double[ m_cCount ];
//	double* uniformV = new double[ m_cCount ];
//
//
//	for ( col=0; col<m_cCols; ++col )
//	{
//		for ( row=0; row<m_cRows; ++row )
//		{
//			At( uniformH, row, col ) = UNIFORM_PLUS_MINUS_ONE;
//			At( uniformV, row, col ) = UNIFORM_PLUS_MINUS_ONE;
//		}
//	}
//
//	// filter with gaussian
//
//	double fConvolvedH, fConvolvedV;
//	double fSampleH, fSampleV;
//	double elasticScale = severityFactor * ::GetPreferences().m_dElasticScaling;
//	int xxx, yyy, xxxDisp, yyyDisp;
//	int iiMid = GAUSSIAN_FIELD_SIZE/2;  // GAUSSIAN_FIELD_SIZE is strictly odd
//
//	for ( col=0; col<m_cCols; ++col )
//	{
//		for ( row=0; row<m_cRows; ++row )
//		{
//			fConvolvedH = 0.0;
//			fConvolvedV = 0.0;
//
//			for ( xxx=0; xxx<GAUSSIAN_FIELD_SIZE; ++xxx )
//			{
//				for ( yyy=0; yyy<GAUSSIAN_FIELD_SIZE; ++yyy )
//				{
//					xxxDisp = col - iiMid + xxx;
//					yyyDisp = row - iiMid + yyy;
//
//					if ( xxxDisp<0 || xxxDisp>=m_cCols || yyyDisp<0 || yyyDisp>=m_cRows )
//					{
//						fSampleH = 0.0;
//						fSampleV = 0.0;
//					}
//					else
//					{
//						fSampleH = At( uniformH, yyyDisp, xxxDisp );
//						fSampleV = At( uniformV, yyyDisp, xxxDisp );
//					}
//
//					fConvolvedH += fSampleH * m_GaussianKernel[ yyy ][ xxx ];
//					fConvolvedV += fSampleV * m_GaussianKernel[ yyy ][ xxx ];
//				}
//			}
//
//			At( m_DispH, row, col ) = elasticScale * fConvolvedH;
//			At( m_DispV, row, col ) = elasticScale * fConvolvedV;
//		}
//	}
//
//	delete[] uniformH;
//	delete[] uniformV;
//
//	// next, the scaling of the image by a random scale factor
//	// Horizontal and vertical directions are scaled independently
//
//	double dSFHoriz = severityFactor * ::GetPreferences().m_dMaxScaling / 100.0 * UNIFORM_PLUS_MINUS_ONE;  // m_dMaxScaling is a percentage
//	double dSFVert = severityFactor * ::GetPreferences().m_dMaxScaling / 100.0 * UNIFORM_PLUS_MINUS_ONE;  // m_dMaxScaling is a percentage
//
//
//	int iMid = m_cRows/2;
//
//	for ( row=0; row<m_cRows; ++row )
//	{
//		for ( col=0; col<m_cCols; ++col )
//		{
//			At( m_DispH, row, col ) += dSFHoriz * ( col-iMid );
//			At( m_DispV, row, col ) -= dSFVert * ( iMid-row );  // negative because of top-down bitmap
//		}
//	}
//
//
//	// finally, apply a rotation
//
//	double angle = severityFactor * ::GetPreferences().m_dMaxRotation * UNIFORM_PLUS_MINUS_ONE;
//	angle = angle * 3.1415926535897932384626433832795 / 180.0;  // convert from degrees to radians
//
//	double cosAngle = cos( angle );
//	double sinAngle = sin( angle );
//
//	for ( row=0; row<m_cRows; ++row )
//	{
//		for ( col=0; col<m_cCols; ++col )
//		{
//			At( m_DispH, row, col ) += ( col-iMid ) * ( cosAngle - 1 ) - ( iMid-row ) * sinAngle;
//			At( m_DispV, row, col ) -= ( iMid-row ) * ( cosAngle - 1 ) + ( col-iMid ) * sinAngle;  // negative because of top-down bitmap
//		}
//	}
//
//
//
//
//	//////////////////////////////////////////////////////////
//	// applies the current distortion map to the input vector
//	
//	// For the mapped array, we assume that 0.0 == background, and 1.0 == full intensity information
//	// This is different from the input vector, in which +1.0 == background (white), and 
//	// -1.0 == information (black), so we must convert one to the other
//	
//	std::vector< std::vector< double > >   mappedVector( m_cRows, std::vector< double >( m_cCols, 0.0 ));
//	
//	double sourceRow, sourceCol;
//	double fracRow, fracCol;
//	double w1, w2, w3, w4;
//	double sourceValue;
//	int row, col;
//	int sRow, sCol, sRowp1, sColp1;
//	BOOL bSkipOutOfBounds;
//	
//	for ( row=0; row<m_cRows; ++row )
//	{
//		for ( col=0; col<m_cCols; ++col )
//		{
//			// the pixel at sourceRow, sourceCol is an "phantom" pixel that doesn't really exist, and
//			// whose value must be manufactured from surrounding real pixels (i.e., since 
//			// sourceRow and sourceCol are floating point, not ints, there's not a real pixel there)
//			// The idea is that if we can calculate the value of this phantom pixel, then its 
//			// displacement will exactly fit into the current pixel at row, col (which are both ints)
//			
//			sourceRow = (double)row - At( m_DispV, row, col );
//			sourceCol = (double)col - At( m_DispH, row, col );
//			
//			// weights for bi-linear interpolation
//			
//			fracRow = sourceRow - (int)sourceRow;
//			fracCol = sourceCol - (int)sourceCol;
//			
//			
//			w1 = ( 1.0 - fracRow ) * ( 1.0 - fracCol );
//			w2 = ( 1.0 - fracRow ) * fracCol;
//			w3 = fracRow * ( 1 - fracCol );
//			w4 = fracRow * fracCol;
//			
//			
//			// limit indexes
//
///*
//			while (sourceRow >= m_cRows ) sourceRow -= m_cRows;
//			while (sourceRow < 0 ) sourceRow += m_cRows;
//			
//			while (sourceCol >= m_cCols ) sourceCol -= m_cCols;
//			while (sourceCol < 0 ) sourceCol += m_cCols;
//*/
//			bSkipOutOfBounds = FALSE;
//
//			if ( (sourceRow + 1.0) >= m_cRows )	bSkipOutOfBounds = TRUE;
//			if ( sourceRow < 0 )				bSkipOutOfBounds = TRUE;
//			
//			if ( (sourceCol + 1.0) >= m_cCols )	bSkipOutOfBounds = TRUE;
//			if ( sourceCol < 0 )				bSkipOutOfBounds = TRUE;
//			
//			if ( bSkipOutOfBounds == FALSE )
//			{
//				// the supporting pixels for the "phantom" source pixel are all within the 
//				// bounds of the character grid.
//				// Manufacture its value by bi-linear interpolation of surrounding pixels
//				
//				sRow = (int)sourceRow;
//				sCol = (int)sourceCol;
//				
//				sRowp1 = sRow + 1;
//				sColp1 = sCol + 1;
//				
//				while (sRowp1 >= m_cRows ) sRowp1 -= m_cRows;
//				while (sRowp1 < 0 ) sRowp1 += m_cRows;
//				
//				while (sColp1 >= m_cCols ) sColp1 -= m_cCols;
//				while (sColp1 < 0 ) sColp1 += m_cCols;
//				
//				// perform bi-linear interpolation
//				
//				sourceValue =	w1 * At( inputVector, sRow  , sCol   ) +
//					w2 * At( inputVector, sRow  , sColp1 ) +
//					w3 * At( inputVector, sRowp1, sCol   ) +
//					w4 * At( inputVector, sRowp1, sColp1 );
//			}
//			else
//			{
//				// At least one supporting pixel for the "phantom" pixel is outside the
//				// bounds of the character grid. Set its value to "background"
//
//				sourceValue = 1.0;  // "background" color in the -1 -> +1 range of inputVector
//			}
//			
//			mappedVector[ row ][ col ] = 0.5 * ( 1.0 - sourceValue );  // conversion to 0->1 range we are using for mappedVector
//			
//		}
//	}
//	
//	// now, invert again while copying back into original vector
//	
//	for ( row=0; row<m_cRows; ++row )
//	{
//		for ( col=0; col<m_cCols; ++col )
//		{
//			At( inputVector, row, col ) = 1.0 - 2.0 * mappedVector[ row ][ col ];
//		}
//	}			
//	
//}



//void MakeSamplesFile(char* SamplesDir, int SamplesNum, char* SamplesFileName)
//{
//	int Dirlen = strlen(SamplesDir);
//	char** FileName = new char*[SamplesNum];
//	for(int i = 0; i < SamplesNum; i++)
//	{
//		FileName[i] = new char[255];
//	}
//
//	int num = GetFileNameFromDir(SamplesDir, FileName, SamplesNum);
//
//	if(num != SamplesNum)
//	{
//		printf("Samples num in Dir:%s is %d, != SamplesNum:%d", SamplesDir, num, SamplesNum);
//		return;
//	}
//
//	//满足开源Demo需要的文件头格式
//	unsigned long  nTest = num;
//	unsigned long  magicImages = 2051;
//	unsigned long  magicLabels = 2049;
//	unsigned long  myrows = 28;
//	nTest = ntohl(nTest);
//	magicImages = ntohl(magicImages);
//	magicLabels = ntohl(magicLabels);
//	myrows = ntohl(myrows);
//
//	CFile my0_5TestImages;
//	CFile my0_5TestLabels;
//
//	//输出文件名字
//	char ImagesName[255];
//	char LabelsName[255];
//
//	sprintf_s(ImagesName, "%sImages%d.idx3-ubyte", SamplesFileName, SamplesNum);
//	sprintf_s(LabelsName, "%sLabels%d.idx1-ubyte", SamplesFileName, SamplesNum);
//
//	//测试样本
//	my0_5TestImages.Open(LPCTSTR(ImagesName), CFile::modeWrite|CFile::modeCreate);
//	my0_5TestImages.Write(&magicImages, sizeof(int));
//	my0_5TestImages.Write(&nTest, sizeof(int));
//	my0_5TestImages.Write(&myrows, sizeof(int));
//	my0_5TestImages.Write(&myrows, sizeof(int));
//	//测试样本标签
//	my0_5TestLabels.Open(LPCTSTR(LabelsName), CFile::modeWrite|CFile::modeCreate);
//	my0_5TestLabels.Write(&magicLabels, sizeof(int));
//	my0_5TestLabels.Write(&nTest, sizeof(int));
//
//	IplImage* pSrc = NULL;
//	IplImage* pNormImg = cvCreateImage(cvSize(28,28), 8, 1);
//
//	int tempLabel;
//	unsigned char tempLabel_char;
//	BYTE * data;
//	unsigned char pArray[ g_cImageSize * g_cImageSize ] = {0};
//	int LabelIdx = Dirlen + 1;
//	
//	int bits[4] = {0};
//	int pow_10[4] = {1,10,100,1000};
//	//循环读入样本
//	for(int i = 0; i < num; i++)
//	{
//		cout<<i<<"\t";
//		pSrc = cvLoadImage(FileName[i], 0);
//
//
//		//样本的标签根据样本的文件名字得到，所以必须在调用前确定相应路径下对应标签所在文件命中的位置
//		//tempLabel = FileName[i][54] - 48;//F:\手写识别项目\字母样本\DW_ZiMu20140127\ZiMu20140127
//		//tempLabel =  FileName[i][29] - 48;//F:\\手写识别项目\\字母样本\\a-f
//		tempLabel_char =  FileName[i][LabelIdx] - 48;
//		
//		//当有label大于9时:
//		//int j = 0;
//		//while(FileName[i][LabelIdx+j] != '_')
//		//{
//		//	bits[j] = FileName[i][LabelIdx+j] - 48;
//		//	j++;
//		//}
//		//tempLabel = 0;
//		//int ii = 0;
//		//while(j)
//		//{
//		//	j--;
//		//	tempLabel += bits[j] * pow_10[ii++];
//		//}
//
//		//对原始图像进行归一化
//		//IplImage* pDst = NULL;
//		//GetCharFromBackground(pSrc, &pDst);
//		//normImg1(pSrc, pNormImg);
//		//cvReleaseImage(&pDst);
//		//cvSaveImage(FileName[i], pNormImg);
//		//continue;
//
//		//cvNamedWindow("After norm", 0);
//		//cvShowImage("After norm", pNormImg);
//		//cvWaitKey();
//
//		data=(BYTE *)pSrc->imageData;
//		for(int ii = 0; ii < g_cImageSize; ii++)
//		{
//			for(int j = 0; j < g_cImageSize; j++)
//			{
//				pArray[ii*28+j] = 255 - data[ii*pSrc->widthStep+j];
//			}
//		}
//
//		my0_5TestImages.Write(pArray, sizeof(pArray));
//		//my0_5TestLabels.Write(&tempLabel, sizeof(int));
//		my0_5TestLabels.Write(&tempLabel_char, 1);
//		cvReleaseImage(&pSrc);
//
//	}
//
//	my0_5TestImages.Close();
//	my0_5TestLabels.Close();
//	cvReleaseImage(&pNormImg);
//}

void normImg1(IplImage * src,IplImage * Des)
{
	int i,j;
	int width=src->width;
	int height=src->height;
	int widestep=src->widthStep;

	//找到最大的高度宽度
	int Ma=(int)(max(width,height));
	int WH=1.2*Ma;

	//申请最大1.44倍的图像
	IplImage * s1=cvCreateImage(cvSize(WH,WH),8,1);

	//初始化为白色的
	memset(s1->imageData,255,s1->height*s1->widthStep);

	int dH=(WH-height)/2;
	int dW=(WH-width)/2;

	//将原图的中心移动到s1中
	for (i=0;i<height;i++)
	{
		for (j=0;j<width;j++)
		{
			s1->imageData[(i+dH)*s1->widthStep+(j+dW)]=src->imageData[i*src->widthStep+j];
		}
	}

	cvResize(s1,Des);

	BYTE * data=(BYTE *)Des->imageData;
	////取反
	//for (i=0;i<28;i++)
	//{
	//	for (j=0;j<28;j++)
	//	{
	//		data[i*Des->widthStep+j]=data[i*Des->widthStep+j];
	//	}
	//}

	cvReleaseImage(&s1);
}


void normImg(IplImage * src,IplImage * Des)
{
	int i,j;
	int width=src->width;
	int height=src->height;
	int widestep=src->widthStep;

	//归一化到20*20的矩形框中
	IplImage * s1=cvCreateImage(cvSize(20,20),8,1);
	cvResize(src,s1);
	BYTE * data=(BYTE *)s1->imageData;

	//计算质心
	int c=0;
	int cx=0;
	int cy=0;
	for (i=0;i<20;i++)
	{
		for (j=0;j<20;j++)
		{
			if(data[i*s1->widthStep+j]<128)
			{
				c++;
				cx+=j;
				cy+=i;
			}
		}
	}

	int Gx=cx/c;
	int Gy=cy/c;

	//根据质心归一化到28*28
	memset(Des->imageData,0,Des->height*Des->widthStep);

	//质心移动到矩形框的中心
	for (i=max(0,14-Gx);i<min(28,14+20-Gx);i++)
	{
		for (j=max(0,14-Gy);j<min(28,14+20-Gy);j++)
		{
			Des->imageData[i*Des->widthStep+j]=data[(i-max(0,14-Gx))*s1->widthStep+j-max(0,14-Gy)];
		}
	}

	cvReleaseImage(&s1);
}

void GetCharFromBackground(IplImage* kSrc, IplImage** pDst)
{
	int left,right,top,bottom;
	int i,j;
	//横向投影
	int * H=new int[kSrc->height];
	memset(H,0,sizeof(int)*kSrc->height);
	int * V=new int[kSrc->width];
	memset(V,0,sizeof(int)*kSrc->width);
	BYTE * Data=(BYTE *)kSrc->imageData;
	for (i=0;i<kSrc->height;i++)
	{
		for (j=0;j<kSrc->width;j++)
		{
			if (Data[i*kSrc->widthStep+j]<128)
			{
				H[i]++;
				V[j]++;
			}
		}
	}

	for (i=0;i<kSrc->height;i++)
	{
		if (H[i]>0)
		{
			top=i;
			break;
		}
	}

	for (i=kSrc->height-1;i>=0;i--)
	{
		if (H[i]>0)
		{
			bottom=i;
			break;
		}
	}

	for (i=0;i<kSrc->width;i++)
	{
		if (V[i]>0)
		{
			left=i;
			break;
		}
	}

	for (i=kSrc->width-1;i>=0;i--)
	{
		if (V[i]>0)
		{
			right=i;
			break;
		}
	}

	CvRect rc;
	rc.x=left;
	rc.y=top;
	rc.width=right-left+1;
	rc.height=bottom-top+1;
	if(*pDst != NULL)
		cvReleaseImage(&(*pDst));
	*pDst=cvCreateImage(cvSize(rc.width,rc.height),kSrc->depth,kSrc->nChannels);
	cvSetImageROI(kSrc,rc);

	cvCopy(kSrc,*pDst);

	cvResetImageROI(kSrc);
	delete []H;
	delete []V;
}
