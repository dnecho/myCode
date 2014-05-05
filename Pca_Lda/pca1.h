#pragma once
void Householder_Tri_Symetry_Diagonal( double a[], int n, double q[],double b[], double c[] );

int Tri_Symmetry_Diagonal_Eigenvector( int n, double b[], double c[], double q[], double eps, int l );

int SymmetricRealMatrix_Eigen( double CovMatrix[], int n,double Eigen[], double EigenVector[] );

int pca1( double * data, int row, int col, double * signals, double * PC, double * V );

int project2PCA( double * newdata, int row, int col,double * PC, double * newsignals, double * ShiftValue );