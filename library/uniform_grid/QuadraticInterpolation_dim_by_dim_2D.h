#ifndef QUADRATIC_INTERPOLATION_2D_H
#define QUADRATIC_INTERPOLATION_2D_H

#include "../ARRAY2D.h"

struct QuadraticInterpolation_dim_by_dim_2D
{
	const ARRAY2D<double>* pU;
	//
	QuadraticInterpolation_dim_by_dim_2D( const ARRAY2D<double>& U )
	{pU=&U;}
	//
	double interpolate( double x, double y )
	{
		const ARRAY2D<double>& U=*pU;
		int isize = U.isize;
		int jsize = U.jsize;
		// select i,j such that (x,y)\in[i,i+1]x[j,j+1]
		// [i-1,i+2]\in[0,isize-1]
		// [j-1,j+2]\in[0,jsize-1]
		int i = (int)(x); if(i<1)i=1; if(i>=isize-2)i=isize-3;
		int j = (int)(y); if(j<1)j=1; if(j>=jsize-2)j=jsize-3;
		// coordinates w.r.t. (i,j)
		x-=i;
		y-=j;
		double xx=x*x,xxx=xx*x;
		double yy=y*y,yyy=yy*y;
		//
		double phi_im1 =  -x/3. +xx/2.-xxx/6., phi_jm1 =  -y/3.+yy/2.-yyy/6.;
		double phi_i   = 1-x/2. -xx   +xxx/2., phi_j   = 1-y/2.-yy   +yyy/2.;
		double phi_ip1 =   x    +xx/2.-xxx/2., phi_jp1 =   y   +yy/2.-yyy/2.;
		double phi_ip2 =  -x/6.       +xxx/6., phi_jp2 =  -y/6.      +yyy/6.;
		//
		return phi_im1*phi_jm1*U[i-1][j-1]
			  +phi_im1*phi_j  *U[i-1][j  ]
			  +phi_im1*phi_jp1*U[i-1][j+1]
			  +phi_im1*phi_jp2*U[i-1][j+2]
			  +phi_i  *phi_jm1*U[i  ][j-1]
			  +phi_i  *phi_j  *U[i  ][j  ]
			  +phi_i  *phi_jp1*U[i  ][j+1]
			  +phi_i  *phi_jp2*U[i  ][j+2]
			  +phi_ip1*phi_jm1*U[i+1][j-1]
			  +phi_ip1*phi_j  *U[i+1][j  ]
			  +phi_ip1*phi_jp1*U[i+1][j+1]
			  +phi_ip1*phi_jp2*U[i+1][j+2]
			  +phi_ip2*phi_jm1*U[i+2][j-1]
			  +phi_ip2*phi_j  *U[i+2][j  ]
			  +phi_ip2*phi_jp1*U[i+2][j+1]
			  +phi_ip2*phi_jp2*U[i+2][j+2];
	}
	//
	double diff_x( double x, double y, double z )
	{
		const ARRAY2D<double>& U=*pU;
		int isize = U.isize;
		int jsize = U.jsize;
		// select i,j such that (x,y)\in[i,i+1]x[j,j+1]
		// [i-1,i+2]\in[0,isize-1]
		// [j-1,j+2]\in[0,jsize-1]
		int i = (int)(x); if(i<1)i=1; if(i>=isize-2)i=isize-3;
		int j = (int)(y); if(j<1)j=1; if(j>=jsize-2)j=jsize-3;
		// coordinates w.r.t. (i,j)
		x-=i;
		y-=j;
		double xx=x*x,xxx=xx*x;
		double yy=y*y,yyy=yy*y;
		//
		double phi_im1 = -1/3.+ x   - xx/2.   , phi_jm1 =  -y/3.+yy/2.-yyy/6.;
		double phi_i   = -1/2.- x*2.+ xx/2.*3., phi_j   = 1-y/2.-yy   +yyy/2.;
		double phi_ip1 =  1   + x   - xx/2.*3., phi_jp1 =   y   +yy/2.-yyy/2.;
		double phi_ip2 = -1/6.      + xx/2.   , phi_jp2 =  -y/6.      +yyy/6.;
		//
		return phi_im1*phi_jm1*U[i-1][j-1]
			  +phi_im1*phi_j  *U[i-1][j  ]
			  +phi_im1*phi_jp1*U[i-1][j+1]
			  +phi_im1*phi_jp2*U[i-1][j+2]
			  +phi_i  *phi_jm1*U[i  ][j-1]
			  +phi_i  *phi_j  *U[i  ][j  ]
			  +phi_i  *phi_jp1*U[i  ][j+1]
			  +phi_i  *phi_jp2*U[i  ][j+2]
			  +phi_ip1*phi_jm1*U[i+1][j-1]
			  +phi_ip1*phi_j  *U[i+1][j  ]
			  +phi_ip1*phi_jp1*U[i+1][j+1]
			  +phi_ip1*phi_jp2*U[i+1][j+2]
			  +phi_ip2*phi_jm1*U[i+2][j-1]
			  +phi_ip2*phi_j  *U[i+2][j  ]
			  +phi_ip2*phi_jp1*U[i+2][j+1]
			  +phi_ip2*phi_jp2*U[i+2][j+2];
	}
	//
	//
	double diff_y( double x, double y )
	{
		const ARRAY2D<double>& U=*pU;
		int isize = U.isize;
		int jsize = U.jsize;
		// select i,j,k such that (x,y,z)\in[i,i+1]x[j,j+1]
		// [i-1,i+2]\in[0,isize-1]
		// [j-1,j+2]\in[0,jsize-1]
		int i = (int)(x); if(i<1)i=1; if(i>=isize-2)i=isize-3;
		int j = (int)(y); if(j<1)j=1; if(j>=jsize-2)j=jsize-3;
		// coordinates w.r.t. (i,j)
		x-=i;
		y-=j;
		double xx=x*x,xxx=xx*x;
		double yy=y*y,yyy=yy*y;
		//
		double phi_im1 =   -x/3. +xx/2.-xxx/6., phi_jm1 = -1/3.+ y   - yy/2.   ;
		double phi_i   =  1-x/2. -xx   +xxx/2., phi_j   = -1/2.- y*2.+ yy/2.*3.;
		double phi_ip1 =    x    +xx/2.-xxx/2., phi_jp1 =  1   + y   - yy/2.*3.;
		double phi_ip2 =   -x/6.       +xxx/6., phi_jp2 = -1/6.      + yy/2.   ;
		//
		return phi_im1*phi_jm1*U[i-1][j-1]
			  +phi_im1*phi_j  *U[i-1][j  ]
			  +phi_im1*phi_jp1*U[i-1][j+1]
			  +phi_im1*phi_jp2*U[i-1][j+2]
			  +phi_i  *phi_jm1*U[i  ][j-1]
			  +phi_i  *phi_j  *U[i  ][j  ]
			  +phi_i  *phi_jp1*U[i  ][j+1]
			  +phi_i  *phi_jp2*U[i  ][j+2]
			  +phi_ip1*phi_jm1*U[i+1][j-1]
			  +phi_ip1*phi_j  *U[i+1][j  ]
			  +phi_ip1*phi_jp1*U[i+1][j+1]
			  +phi_ip1*phi_jp2*U[i+1][j+2]
			  +phi_ip2*phi_jm1*U[i+2][j-1]
			  +phi_ip2*phi_j  *U[i+2][j  ]
			  +phi_ip2*phi_jp1*U[i+2][j+1]
			  +phi_ip2*phi_jp2*U[i+2][j+2];
	}
};

#endif
