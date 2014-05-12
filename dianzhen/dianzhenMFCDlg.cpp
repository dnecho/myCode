// dianzhenMFCDlg.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "dianzhenMFC.h"
#include "dianzhenMFCDlg.h"

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


// CdianzhenMFCDlg �Ի���




CdianzhenMFCDlg::CdianzhenMFCDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CdianzhenMFCDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CdianzhenMFCDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CdianzhenMFCDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BUTTON1, &CdianzhenMFCDlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_BUTTON2, &CdianzhenMFCDlg::OnBnClickedButton2)
END_MESSAGE_MAP()


// CdianzhenMFCDlg ��Ϣ�������

BOOL CdianzhenMFCDlg::OnInitDialog()
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
	pSrc3C = NULL;
	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
}

void CdianzhenMFCDlg::OnSysCommand(UINT nID, LPARAM lParam)
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

void CdianzhenMFCDlg::OnPaint()
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
HCURSOR CdianzhenMFCDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void CdianzhenMFCDlg::OnBnClickedButton1()
{
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
		}

		pSrc = cvLoadImage(FileName.GetBuffer(0),0);
		pSrc3C = cvLoadImage(FileName.GetBuffer(0),1);

		//int dwidth=500;
		//int dheight=dwidth*pSrc->height/pSrc->width;

		//psrc=cvCreateImage(cvSize(dwidth,dheight),src->depth,src->nChannels);
		//cvResize(src,psrc);

		//cvNamedWindow("ԭͼ",0);
		//cvShowImage("ԭͼ",pSrc3C);

		//cvWaitKey();

	}
}

void CdianzhenMFCDlg::OnBnClickedButton2()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	if(pSrc==NULL||pSrc3C==NULL)
	{
		AfxMessageBox("���ͼƬ");
		return;
	}
	
	long t1=GetTickCount();
	
	Detect(pSrc, pSrc3C);
	long t2=GetTickCount();

	CString ProjectTime;
	ProjectTime.Format("%d", t2-t1);
	AfxMessageBox(ProjectTime);
		
	cvNamedWindow("ʶ����", 0);
	cvShowImage("ʶ����", pSrc3C);
	cvWaitKey();

}
