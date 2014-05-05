#pragma once

void ppp( double a[], double e[], double s[],double v[], int m, int n );

void sss( double fg[2], double cs[2] );

int svd( double a[], int m, int n, double u[], double v[], double eps, int ka );

int pca2( double * data, int row, int col, double * signals, double * PC, double * V );