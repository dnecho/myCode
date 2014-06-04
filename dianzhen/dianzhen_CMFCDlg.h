// dianzhen_CMFCDlg.h : 头文件
//

#pragma once


// Cdianzhen_CMFCDlg 对话框
class Cdianzhen_CMFCDlg : public CDialog
{
// 构造
public:
	Cdianzhen_CMFCDlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_DIANZHEN_CMFC_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持

	IplImage* pSrc;
	CString PicName;
	CString Outcome;

// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButton1();
	afx_msg void OnBnClickedButton2();
};
