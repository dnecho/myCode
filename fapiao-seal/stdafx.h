// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>



// TODO: 在此处引用程序需要的其他头文件
#include "opencv2/opencv.hpp"
using namespace std;
using namespace cv;
#pragma comment(lib,"cvLib/opencv_video231.lib")
#pragma comment(lib,"cvLib/opencv_video231d.lib")
#pragma comment(lib,"cvLib/opencv_highgui231.lib")
#pragma comment(lib,"cvLib/opencv_highgui231d.lib")
#pragma comment(lib,"cvLib/opencv_core231.lib")
#pragma comment(lib,"cvLib/opencv_core231d.lib")
#pragma comment(lib,"cvLib/opencv_imgproc231.lib")
#pragma comment(lib,"cvLib/opencv_imgproc231d.lib")
#pragma comment(lib,"cvLib/opencv_objdetect231.lib")
#pragma comment(lib,"cvLib/opencv_objdetect231d.lib")


#pragma comment(lib,"cvLib/opencv_ml231.lib")
#pragma comment(lib,"cvLib/opencv_ml231d.lib")
#include "region.h"
#include "rotate.h"