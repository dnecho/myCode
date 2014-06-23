// dianzhen_CMFCDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "dianzhen_CMFC.h"
#include "dianzhen_CMFCDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// 对话框数据
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// Cdianzhen_CMFCDlg 对话框




Cdianzhen_CMFCDlg::Cdianzhen_CMFCDlg(CWnd* pParent /*=NULL*/)
	: CDialog(Cdianzhen_CMFCDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void Cdianzhen_CMFCDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(Cdianzhen_CMFCDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BUTTON1, &Cdianzhen_CMFCDlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_BUTTON2, &Cdianzhen_CMFCDlg::OnBnClickedButton2)
END_MESSAGE_MAP()


// Cdianzhen_CMFCDlg 消息处理程序

BOOL Cdianzhen_CMFCDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: 在此添加额外的初始化代码
	pSrc = NULL;
	PicScale = 1;
	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void Cdianzhen_CMFCDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void Cdianzhen_CMFCDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR Cdianzhen_CMFCDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void Cdianzhen_CMFCDlg::OnBnClickedButton1()
{
	
	// TODO: 在此添加控件通知处理程序代码
	// TODO: 在此添加控件通知处理程序代码
	CFileDialog dlg(TRUE,//TRUE是创建打开文件对话框，FALSE则创建的是保存文件对话框 
		".jpg",//默认的打开文件的类型 
		NULL,//默认打开的文件名 
		OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,//打开只读文件 
		"图片文件|*.jpg;*.png;*.bmp|所有文件 (*.*)|*.*||");//所有可以打开的文件类型 

	if(dlg.DoModal()==IDOK)
	{
		CString FileName=dlg.GetPathName();//全路径名

		if (!pSrc)
		{
			cvReleaseImage(&pSrc);
			pSrc=NULL;
		}

		IplImage* pSrc1 = cvLoadImage(FileName.GetBuffer(0),0);

		PicName=FileName;
		int dwidth=800;
		
		PicScale = dwidth*1.0/pSrc1->width;
		int dheight=(int)(pSrc1->height*PicScale);

		pSrc=cvCreateImage(cvSize(dwidth,dheight),pSrc1->depth,pSrc1->nChannels);
		cvResize(pSrc1,pSrc);
		cvReleaseImage(&pSrc1);
		//cvNamedWindow("原图",0);
		//cvShowImage("原图",pSrc3C);

		//cvWaitKey();

	}
}

void Cdianzhen_CMFCDlg::OnBnClickedButton2()
{
	// TODO: 在此添加控件通知处理程序代码
	Outcome="";
	if(pSrc == NULL)
	{
		AfxMessageBox("请先打开图片");
		return;
	}
	Axis center={0};
	int pixelsize=0;

	unsigned char* pSrcData = (unsigned char*)pSrc->imageData;

	//cvSmooth(pSrc, pSrc, CV_GAUSSIAN, 3);

	unsigned char* pImgData = (unsigned char*)malloc(sizeof(unsigned char)*(pSrc->width*pSrc->height));
	for(int i = 0; i < pSrc->height; i++)
		for(int j = 0; j < pSrc->width; j++)
			pImgData[j+i*pSrc->width] = pSrcData[j+i*pSrc->widthStep];			

	int IsSuccess = process(pImgData, pSrc->width, pSrc->height, &center, &pixelsize, Outcome);
	
	if(IsSuccess==0)
	{
		CString tmp;
		tmp.Format("处理失败，请记录该图片，图片名字为：%s",PicName);
		Outcome += tmp;
		AfxMessageBox(Outcome);
	}else
	{
		CString tmp;
		tmp.Format("求平均值final\tx(单位0.01mm):%d，y(单位0.01mm):%d，PixelSize(单位0.001mm):%d\n", center.x, center.y, (int)(pixelsize*PicScale+0.5));
		Outcome += tmp;
		AfxMessageBox(Outcome);
	}
}
