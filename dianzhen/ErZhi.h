#pragma once

//#include "opencv_IpiImage.h"

void adaptiveThreshold(IplImage * GrayImg, IplImage * BiImg);

void SML_OTSU(IplImage * Src,IplImage * Des);

void SML_OTSU_Bersen(IplImage * Src,IplImage * Des);