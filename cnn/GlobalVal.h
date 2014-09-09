#pragma once

typedef unsigned int UINT;
typedef const char * LPCTSTR;
typedef unsigned char BYTE;
typedef unsigned long DWORD;

const unsigned int g_cImageSize = 28;

//º¯ÊýCMNistDoc::CMNistDoc()ÖÐ
const int m_cCols = 29;
const int m_cRows = 29;
const int m_cCount = m_cCols * m_cRows;

#define GAUSSIAN_FIELD_SIZE ( 21 )  // strictly odd number
#define ASSERT assert

#include <stdio.h>
#include <tchar.h>
#include <stdio.h>
#include <tchar.h>
#include <fstream>
#include <iostream>
#include <io.h>
#include <stdlib.h>
#include <assert.h>

#include "opencv2/opencv.hpp"
using namespace std;
using namespace cv;
#pragma comment(lib,"cvLib/opencv_highgui231.lib")
#pragma comment(lib,"cvLib/opencv_highgui231d.lib")
#pragma comment(lib,"cvLib/opencv_core231.lib")
#pragma comment(lib,"cvLib/opencv_core231d.lib")
#pragma comment(lib,"cvLib/opencv_imgproc231.lib")
#pragma comment(lib,"cvLib/opencv_imgproc231d.lib")

#include "NeuralNetwork.h"
#include "myFun.h"
