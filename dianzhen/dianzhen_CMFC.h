// dianzhen_CMFC.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// Cdianzhen_CMFCApp:
// �йش����ʵ�֣������ dianzhen_CMFC.cpp
//

class Cdianzhen_CMFCApp : public CWinApp
{
public:
	Cdianzhen_CMFCApp();

// ��д
	public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern Cdianzhen_CMFCApp theApp;