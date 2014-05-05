/*
   PCA program (using SVD)
   Written by Y. Bin Mao
   Video and Image Processing and Analysis Group (VIPAG)
   School of Automation, NJUST
   Jan. 8, 2008

   All rights reserved. (c)

   Here is a matlab description of the algorithm
   % PCA2: Perform PCA using SVD.
   % data     --- MxN matrix of input data ( M dimensions, N trials )
   % signals  --- MxN matrix of projected data 
   % PC       --- each column is a PC
   % V        --- Mx1 matrix of variances
   %
   function [signals, PC, V] = pca2( data )

   [M, N] = size( data );

   % subtract off the mean for each dimension
   mn = mean( data, 2 );
   data = data - repmat( mn, 1, N );

   % construct the matrix Y
   Y = data' / sqrt(N-1);

   % SVD does it all
   [u, S, PC] = svd( Y );

   % calculate the variances
   S = diag( S );
   V = S .* S;

   % project the original data
   signals = PC' * data;
*/
#include "stdafx.h"
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

void ppp( double a[], double e[], double s[],double v[], int m, int n )
{ 
	int i,j,p,q;
	double d;
	
	if ( m >= n ) i = n;
	else i = m;
	for ( j = 1; j <= i-1; j++ )
	{ 
		a[(j-1)*n+j-1] = s[j-1];
		a[(j-1)*n+j] = e[j-1];
	}
	a[(i-1)*n+i-1] = s[i-1];
	if ( m < n ) a[(i-1)*n+i] = e[i-1];
	for ( i = 1; i <= n-1; i++ )
		for ( j = i+1; j <= n; j++ )
		{ 
			p = (i-1)*n+j-1; q = (j-1)*n+i-1;
			d = v[p]; v[p] = v[q]; v[q] = d;
		}
		return;
}

void sss( double fg[2], double cs[2] )
{
	double r, d;
    if ( ( fabs(fg[0]) + fabs(fg[1] ) ) == 0.0 )
	{ 
		cs[0] = 1.0; cs[1] = 0.0; d = 0.0;
	}
    else
	{ 
		d = sqrt( fg[0]*fg[0]+fg[1]*fg[1] );
		if ( fabs( fg[0] ) > fabs( fg[1] ) )
		{
			d = fabs(d);
			if ( fg[0] < 0.0 ) d = -d;
		}
        if ( fabs( fg[1] ) >= fabs( fg[0] ) )
        { 
			d = fabs(d);
            if ( fg[1] < 0.0 ) d = -d;
		}
        cs[0] = fg[0]/d; cs[1] = fg[1]/d;
    }
    r = 1.0;
    if ( fabs( fg[0] ) > fabs( fg[1] ) ) 
		r = cs[1];
    else
		if ( cs[0] != 0.0 ) r = 1.0/cs[0];
		fg[0] = d; fg[1] = r;
		
		return;
}

/*
   一般实矩阵奇异值分解
   徐士良. 常用算法程序集（C语言描述），第3版. 清华大学出版社. 2004

   double a[m][n] --- 存放mxn维实矩阵A。返回时其对角线给出奇异值（以非递增次序排列）
                      其余元素均为0
   int m, int n   --- 实矩阵A的行数和列数
   double u[m][m] --- 返回左奇异向量U
   double v[n][n] --- 返回右奇异向量V'
   double eps     --- 给定的精度要求
   int ka         --- 其值为max(m, n) + 1
   返回值：
           若返回标志值小于0，则表示程序工作失败；
		   若返回标志值大于0，则表示正常返回   
*/
int svd( double a[], int m, int n, double u[], double v[], double eps, int ka )
{ 
	int i, j, k, l, it, ll, kk, ix, iy, mm, nn, iz, m1, ks;
    double d, dd, t, sm, sm1, em1, sk, ek, b, c, shh, fg[2], cs[2];
    double * s, * e, * w;

    s = ( double * )malloc( ka * sizeof(double) );
    e = ( double * )malloc( ka * sizeof(double) );
    w = ( double * )malloc( ka * sizeof(double) );
    it = 60; k = n;
	if ( m-1 < n ) k = m-1;
	l = m;
	if ( n-2 < m ) l = n-2;
	if ( l < 0 ) l = 0;
	ll = k;
	if ( l > k ) ll = l;
	if ( ll >= 1 )
	{ 
		for ( kk = 1; kk <= ll; kk++ )
		{ 
			if ( kk <= k )
			{ 
				d = 0.0;
				for ( i = kk; i <= m; i++ )
				{ 
					ix = (i-1)*n+kk-1; d = d+a[ix]*a[ix]; 
				}
                s[kk-1] = sqrt(d);
                if ( s[kk-1] != 0.0 )
				{ 
					ix = (kk-1)*n+kk-1;
                    if ( a[ix] != 0.0 )
					{ 
						s[kk-1] = fabs( s[kk-1] );
                        if ( a[ix] < 0.0 ) s[kk-1] = -s[kk-1];
                    }
                    for ( i = kk; i <= m; i++ )
					{  
						iy = (i-1)*n+kk-1;
                        a[iy] = a[iy]/s[kk-1];
                    }
                    a[ix] = 1.0+a[ix];
                }
                s[kk-1] = -s[kk-1];
            }
            if ( n >= kk+1 )
            { 
				for ( j = kk+1; j <= n; j++ )
                { 
					if (( kk <= k ) && ( s[kk-1] != 0.0 ) )
                    { 
						d = 0.0;
                        for ( i = kk; i <= m; i++ )
                        { 
							ix = (i-1)*n+kk-1;
                            iy = (i-1)*n+j-1;
                            d = d+a[ix]*a[iy];
                        }
                        d = -d/a[(kk-1)*n+kk-1];
                        for ( i = kk; i <= m; i++ )
                        { 
							ix = (i-1)*n+j-1;
                            iy = (i-1)*n+kk-1;
                            a[ix] = a[ix]+d*a[iy];
                        }
                    }
                    e[j-1] = a[(kk-1)*n+j-1];
                }
            }
            if ( kk <= k )
            { 
				for ( i = kk; i <= m; i++ )
                { 
					ix = (i-1)*m+kk-1; 
					iy = (i-1)*n+kk-1;
                    u[ix] = a[iy];
                }
            }
            if ( kk <= l )
            { 
				d = 0.0;
                for ( i = kk+1; i <= n; i++ )
					d = d + e[i-1]*e[i-1];
                e[kk-1] = sqrt(d);
                if ( e[kk-1] != 0.0 )
                { 
					if ( e[kk] != 0.0 )
                    { 
						e[kk-1] = fabs(e[kk-1]);
                        if ( e[kk] < 0.0 ) e[kk-1] = -e[kk-1];
                    }
                    for ( i = kk+1; i <= n; i++ )
						e[i-1] = e[i-1]/e[kk-1];
                    e[kk] = 1.0+e[kk];
                }
                e[kk-1] = -e[kk-1];
                if (( kk+1 <= m ) && ( e[kk-1] != 0.0 ) )
                { 
					for ( i = kk+1; i <= m; i++ ) w[i-1] = 0.0;
                    for ( j = kk+1; j <= n; j++)
						for ( i = kk+1; i <= m; i++ )
							w[i-1] = w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for ( j = kk+1; j <= n; j++ )
                        for ( i = kk+1; i <= m; i++ )
                        { 
							ix = (i-1)*n+j-1;
                            a[ix] = a[ix]-w[i-1]*e[j-1]/e[kk];
                        }
                }
                for ( i = kk+1; i <= n; i++ )
					v[(i-1)*n+kk-1] = e[i-1];
            }
        }
	}
    mm = n;
    if ( m+1 < n ) mm = m+1;
    if ( k < n ) s[k] = a[k*n+k];
    if ( m < mm ) s[mm-1] = 0.0;
    if ( l+1 < mm) e[l] = a[l*n+mm-1];
    e[mm-1] = 0.0;
    nn = m;
    if ( m > n ) nn = n;
    if ( nn >= k+1 )
    { 
		for ( j = k+1; j <= nn; j++ )
        { 
			for ( i = 1; i <= m; i++ )
				u[(i-1)*m+j-1] = 0.0;
            u[(j-1)*m+j-1] = 1.0;
        }
    }
    if ( k >= 1)
    { 
		for (ll = 1; ll <= k; ll++ )
        { 
			kk = k-ll+1; iz = (kk-1)*m+kk-1;
            if ( s[kk-1] != 0.0 )
            {
				if ( nn >= kk+1 )
                  for ( j = kk+1; j <= nn; j++ )
                  { 
					  d = 0.0;
                      for ( i = kk; i <= m; i++ )
                      {   
						  ix = (i-1)*m+kk-1;
                          iy = (i-1)*m+j-1;
                          d = d+u[ix]*u[iy]/u[iz];
                      }
                      d = -d;
                      for ( i = kk; i <= m; i++ )
                      { 
						  ix = (i-1)*m+j-1;
                          iy = (i-1)*m+kk-1;
                          u[ix] = u[ix]+d*u[iy];
                      }
                  }
                  for ( i = kk; i <= m; i++ )
                  { 
					  ix = (i-1)*m+kk-1; u[ix] = -u[ix];
				  }
                  u[iz] = 1.0+u[iz];
                  if ( kk-1 >= 1 )
					  for ( i = 1; i <= kk-1; i++ )
						  u[(i-1)*m+kk-1] = 0.0;
            }
            else
            { 
				for ( i = 1; i <= m; i++ )
					u[(i-1)*m+kk-1] = 0.0;
                u[(kk-1)*m+kk-1] = 1.0;
            }
		}
    }
    for ( ll = 1; ll <= n; ll++ )
	{ 
		kk = n-ll+1; iz = kk*n+kk-1;
        if ( (kk<=l) && (e[kk-1] != 0.0) )
        { 
			for ( j = kk+1; j <= n; j++ )
            { 
				d = 0.0;
                for ( i = kk+1; i <= n; i++ )
                { 
					ix = (i-1)*n+kk-1; iy = (i-1)*n+j-1;
                    d = d+v[ix]*v[iy]/v[iz];
                }
                d = -d;
                for ( i = kk+1; i <= n; i++ )
                { 
					ix = (i-1)*n+j-1; iy = (i-1)*n+kk-1;
                    v[ix] = v[ix]+d*v[iy];
                }
            }
        }
        for ( i = 1; i <= n; i++ )
			v[(i-1)*n+kk-1] = 0.0;
        v[iz-n] = 1.0;
    }
    for ( i = 1; i <= m; i++ )
    for ( j = 1; j <= n; j++ )
		a[(i-1)*n+j-1] = 0.0;
    m1 = mm; it = 60;
    while ( 1==1 )
    {
		if ( mm == 0 )
        { 
			ppp( a, e, s, v, m, n );
            free( s ); 
			free( e ); 
			free( w ); 
			return( 1 );
        }
        if ( it == 0 )
        { 
			ppp( a, e, s, v, m, n );
            free( s ); 
			free( e ); 
			free( w ); 
			return( -1 );
        }
        kk = mm-1;
		while ( (kk != 0) && (fabs(e[kk-1]) != 0.0) )
		{ 
			d = fabs(s[kk-1])+fabs(s[kk]);
			dd = fabs(e[kk-1]);
			if ( dd > eps*d ) 
				kk = kk-1;
			else 
				e[kk-1] = 0.0;
		}
		if ( kk == mm-1 )
		{ 
			kk = kk+1;
			if ( s[kk-1] < 0.0 )
			{ 
				s[kk-1] = -s[kk-1];
				for ( i = 1; i <= n; i++ )
				{ 
					ix = (i-1)*n+kk-1; 
					v[ix] = -v[ix];
				}
			}
			while ( (kk!=m1) && (s[kk-1] < s[kk]) )
			{ 
				d = s[kk-1]; 
				s[kk-1] = s[kk]; 
				s[kk] = d;
				if ( kk < n )
					for (i=1; i<=n; i++)
					{ 
						ix = (i-1)*n+kk-1; iy = (i-1)*n+kk;
						d = v[ix]; v[ix] = v[iy]; v[iy] = d;
					}
					if ( kk < m )
						for ( i = 1; i <= m; i++ )
						{ 
							ix = (i-1)*m+kk-1; iy = (i-1)*m+kk;
							d = u[ix]; u[ix] = u[iy]; u[iy] = d;
						}
					kk = kk+1;
			}
			it = 60;
			mm = mm-1;
		}
		else
		{ 
			ks = mm;
			while ( (ks>kk) && (fabs(s[ks-1]) != 0.0 ) )
			{ 
				d = 0.0;
				if ( ks != mm ) d = d+fabs(e[ks-1]);
				if ( ks != kk+1 ) d = d+fabs(e[ks-2]);
				dd = fabs(s[ks-1]);
				if ( dd > eps*d ) ks = ks-1;
				else s[ks-1] = 0.0;
			}
			if ( ks == kk )
			{ 
				kk = kk+1;
				d = fabs(s[mm-1]);
				t = fabs(s[mm-2]);
				if ( t > d ) d = t;
				t = fabs( e[mm-2] );
				if ( t > d ) d = t;
				t = fabs( s[kk-1] );
				if ( t > d ) d = t;
				t = fabs( e[kk-1] );
				if ( t > d ) d = t;
				sm = s[mm-1]/d; sm1 = s[mm-2]/d;
				em1 = e[mm-2]/d;
				sk = s[kk-1]/d; ek = e[kk-1]/d;
				b = ((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
				c = sm*em1; c = c*c; shh = 0.0;
				if ( (b!=0.0) || (c!=0.0) )
				{ 
					shh = sqrt( b*b+c );
					if ( b < 0.0 ) shh = -shh;
					shh = c/(b+shh);
				}
				fg[0] = (sk+sm)*(sk-sm)-shh;
				fg[1] = sk*ek;
				for ( i = kk; i <= mm-1; i++ )
				{ 
					sss( fg, cs );
					if ( i != kk ) e[i-2] = fg[0];
					fg[0] = cs[0]*s[i-1]+cs[1]*e[i-1];
					e[i-1] = cs[0]*e[i-1]-cs[1]*s[i-1];
					fg[1] = cs[1]*s[i];
					s[i] = cs[0]*s[i];
					if ( (cs[0] != 1.0) || (cs[1] != 0.0) )
						for ( j = 1; j <= n; j++ )
						{ 
							ix = (j-1)*n+i-1;
							iy = (j-1)*n+i;
							d = cs[0]*v[ix]+cs[1]*v[iy];
							v[iy] = -cs[1]*v[ix]+cs[0]*v[iy];
							v[ix] = d;
						}
						sss( fg, cs );
						s[i-1] = fg[0];
						fg[0] = cs[0]*e[i-1]+cs[1]*s[i];
						s[i] = -cs[1]*e[i-1]+cs[0]*s[i];
						fg[1] = cs[1]*e[i];
						e[i] = cs[0]*e[i];
						if ( i < m )
							if ( ( cs[0] != 1.0 ) || ( cs[1] != 0.0 ) )
								for ( j = 1; j <= m; j++ )
								{ 
									ix = (j-1)*m+i-1;
									iy = (j-1)*m+i;
									d = cs[0]*u[ix]+cs[1]*u[iy];
									u[iy] = -cs[1]*u[ix]+cs[0]*u[iy];
									u[ix] = d;
								}
				}
				e[mm-2] = fg[0];
				it = it-1;
			}
			else
			{ 
				if ( ks == mm )
				{ 
					kk = kk+1;
					fg[1] = e[mm-2]; e[mm-2] = 0.0;
	                for ( ll = kk; ll <= mm-1; ll++ )
		            { 
						i = mm+kk-ll-1;
				        fg[0] = s[i-1];
					    sss( fg, cs );
						s[i-1] = fg[0];
					    if ( i != kk )
						{ 
							fg[1] = -cs[1]*e[i-2];
							e[i-2] = cs[0]*e[i-2];
						}
						if ( (cs[0]!=1.0)||(cs[1]!=0.0) )
							for ( j = 1; j <= n; j++ )
							{ 
								ix = (j-1)*n+i-1;
								iy = (j-1)*n+mm-1;
								d = cs[0]*v[ix]+cs[1]*v[iy];
								v[iy] = -cs[1]*v[ix]+cs[0]*v[iy];
								v[ix] = d;
							}
					}
				}
				else
				{ 
					kk = ks+1;
					fg[1] = e[kk-2];
					e[kk-2] = 0.0;
					for ( i = kk; i <= mm; i++ )
					{ 
						fg[0] = s[i-1];
						sss( fg, cs );
						s[i-1] = fg[0];
						fg[1] = -cs[1]*e[i-1];
						e[i-1] = cs[0]*e[i-1];
						if ( (cs[0]!=1.0)||(cs[1]!=0.0) )
							for (j=1; j<=m; j++)
							{ 
								ix = (j-1)*m+i-1;
								iy = (j-1)*m+kk-2;
								d = cs[0]*u[ix]+cs[1]*u[iy];
								u[iy] = -cs[1]*u[ix]+cs[0]*u[iy];
								u[ix] = d;
							}
					}
				}
			}
		}
    }

    return( 1 );
}

# define EPS  0.000001

/*
  PCA2: Perform PCA using SVD.
  data     --- MxN matrix of input data ( M dimensions, N trials )
  signals  --- MxN matrix of projected data 
  PC       --- each column is a PC
  V        --- Mx1 matrix of variances
  row = M dimensions, col = N trials 
*/
int pca2( double * data, int row, int col, double * signals, double * PC, double * V )
{
	int x, y, i, ka, rvalue;
	double * mean;
	double * u;
    double * d;
    double * v;
	double * sd;
	double * data1;
	double sqrt_col_1;
	int col1, row1;

	if ( row <= 1 || col <= 1 ) return( -1 );

	/* subtract off the mean for each dimension */
	mean = ( double * )malloc( sizeof( double )*row );
	data1 = ( double * ) malloc( sizeof( double )*row*col );
	for ( y = 0; y < row; y++ ) /* calculate mean */
	{
		mean[y] = 0;
		for ( x = 0; x < col; x++ )
			mean[y] += data[y*col+x];
	}
	for ( y = 0; y < row; y++ ) mean[y] = mean[y]/col;
	for ( y = 0; y < row; y++ ) /* subtract mean */
		for ( x = 0; x < col; x++ )
		{
			data1[y*col+x] = data[y*col+x]-mean[y];
		}
	free( mean );

	/* construct the matrix Y: Y = data' / sqrt(col-1); */

	/* construct the matrix Y: Y = data' / sqrt(col-1); */
	sqrt_col_1 = sqrt(col-1.0);
	row1 = col;
	col1 = row;
	sd = (double *)malloc( sizeof(double)*col1*row1 );
	u = (double *)malloc( sizeof(double)*row1*row1 );
	v = (double *)malloc( sizeof(double)*col1*col1 );
	for ( y = 0; y < row1; y++ )
		for ( x = 0; x < col1; x++ )
		{
			sd[y*col1+x] = data1[x*row1+y]/sqrt_col_1;
		}

	/* SVD does it all: [u, S, PC] = svd( Y ); */
	if ( row1 >= col1 ) 
		ka = row1+1;
	else 
		ka = col1+1;
	rvalue = svd( sd, row1, col1, u, v, EPS, ka ); /* svd decomposition */
	d = (double *)malloc( sizeof(double) * col1 );
	for ( i = 0; i < col1; i++ ) d[i] = 0;
	if ( row1 <= col1 )
	{
		for ( i = 0; i < row1; i++ )
			d[i] = sd[i*col1+i];
	}
	else
	{
		for ( i = 0; i < col1; i++ )
			d[i] = sd[i*col1+i];
	}

	/* calculate the variances: S = diag( S ); V = S .* S; */
	for ( x = 0; x < col1; x++ )
		V[x] = d[x] * d[x];
	/* get PC */
	for ( y = 0; y < col1; y++ )
		for ( x = 0; x < col1; x++ )
		{
			PC[y*col1+x] = v[x*col1+y];
		}

	/* project the original data: signals = PC' * data; */
	for ( y = 0; y < row; y++ )
		for ( x = 0; x < col; x++ )
		{
			signals[y*col+x] = 0;
			for ( i = 0; i < row; i++ )
				signals[y*col+x] += PC[i*row+y] * data1[i*col+x];
		}

	free( u );
	free( v );
	free( sd );
	free( d );
	free( data1 );

	return( 1 );
}