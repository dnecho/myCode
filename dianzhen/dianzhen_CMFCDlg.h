// dianzhen_CMFCDlg.h : ͷ�ļ�
//

#pragma once


// Cdianzhen_CMFCDlg �Ի���
class Cdianzhen_CMFCDlg : public CDialog
{
// ����
public:
	Cdianzhen_CMFCDlg(CWnd* pParent = NULL);	// ��׼���캯��

// �Ի�������
	enum { IDD = IDD_DIANZHEN_CMFC_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV ֧��

	IplImage* pSrc;
	CString PicName;
	CString Outcome;

// ʵ��
protected:
	HICON m_hIcon;

	// ���ɵ���Ϣӳ�亯��
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButton1();
	afx_msg void OnBnClickedButton2();
};
