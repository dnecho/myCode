// dianzhen_CMFCDlg.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "dianzhen_CMFC.h"
#include "dianzhen_CMFCDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// ����Ӧ�ó��򡰹��ڡ��˵���� CAboutDlg �Ի���

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// �Ի�������
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

// ʵ��
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


// Cdianzhen_CMFCDlg �Ի���




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


// Cdianzhen_CMFCDlg ��Ϣ�������

BOOL Cdianzhen_CMFCDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// ��������...���˵�����ӵ�ϵͳ�˵��С�

	// IDM_ABOUTBOX ������ϵͳ���Χ�ڡ�
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

	// ���ô˶Ի����ͼ�ꡣ��Ӧ�ó��������ڲ��ǶԻ���ʱ����ܽ��Զ�
	//  ִ�д˲���
	SetIcon(m_hIcon, TRUE);			// ���ô�ͼ��
	SetIcon(m_hIcon, FALSE);		// ����Сͼ��

	// TODO: �ڴ���Ӷ���ĳ�ʼ������
	pSrc = NULL;
	PicScale = 1;
	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
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

// �����Ի��������С����ť������Ҫ����Ĵ���
//  �����Ƹ�ͼ�ꡣ����ʹ���ĵ�/��ͼģ�͵� MFC Ӧ�ó���
//  �⽫�ɿ���Զ���ɡ�

void Cdianzhen_CMFCDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // ���ڻ��Ƶ��豸������

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// ʹͼ���ڹ����������о���
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// ����ͼ��
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

//���û��϶���С������ʱϵͳ���ô˺���ȡ�ù��
//��ʾ��
HCURSOR Cdianzhen_CMFCDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void Cdianzhen_CMFCDlg::OnBnClickedButton1()
{
	
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	CFileDialog dlg(TRUE,//TRUE�Ǵ������ļ��Ի���FALSE�򴴽����Ǳ����ļ��Ի��� 
		".jpg",//Ĭ�ϵĴ��ļ������� 
		NULL,//Ĭ�ϴ򿪵��ļ��� 
		OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,//��ֻ���ļ� 
		"ͼƬ�ļ�|*.jpg;*.png;*.bmp|�����ļ� (*.*)|*.*||");//���п��Դ򿪵��ļ����� 

	if(dlg.DoModal()==IDOK)
	{
		CString FileName=dlg.GetPathName();//ȫ·����

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
		//cvNamedWindow("ԭͼ",0);
		//cvShowImage("ԭͼ",pSrc3C);

		//cvWaitKey();

	}
}

void Cdianzhen_CMFCDlg::OnBnClickedButton2()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	Outcome="";
	if(pSrc == NULL)
	{
		AfxMessageBox("���ȴ�ͼƬ");
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
		tmp.Format("����ʧ�ܣ����¼��ͼƬ��ͼƬ����Ϊ��%s",PicName);
		Outcome += tmp;
		AfxMessageBox(Outcome);
	}else
	{
		CString tmp;
		tmp.Format("��ƽ��ֵfinal\tx(��λ0.01mm):%d��y(��λ0.01mm):%d��PixelSize(��λ0.001mm):%d\n", center.x, center.y, (int)(pixelsize*PicScale+0.5));
		Outcome += tmp;
		AfxMessageBox(Outcome);
	}
}
