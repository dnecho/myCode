// dianzhenMFC.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CdianzhenMFCApp:
// �йش����ʵ�֣������ dianzhenMFC.cpp
//

class CdianzhenMFCApp : public CWinApp
{
public:
	CdianzhenMFCApp();

// ��д
	public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CdianzhenMFCApp theApp;