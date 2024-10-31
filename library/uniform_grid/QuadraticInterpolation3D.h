#ifndef QUADRATIC_3D_H
#define QUADRATIC_3D_H

#include "../Macros.h"
#include "GRID3D.h"
#include "../ARRAY3D.h"

class QuadraticInterpolation3D : public Function_3D {
protected:
	GRID3D grid;
	const ARRAY3D<double>* pf;
	      ARRAY3D<double> fxx,fyy,fzz;
public:
	QuadraticInterpolation3D( const GRID3D& grid_, const ARRAY3D<double>& f ) {
		grid=grid_; int isize=grid.isize;
		pf=&f;      int jsize=grid.jsize;
		            int ksize=grid.ksize;
		// 
		fxx.resize(isize,jsize,ksize);
		fyy.resize(isize,jsize,ksize);
		fzz.resize(isize,jsize,ksize);
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++)
		for(int k=0;k<ksize;k++){
			double f000 = f[i][j][k];
			double fx_iph = (i==isize-1) ? 0 : f[i+1][j][k]-f000, fx_imh = (i==0) ? 0 : f000-f[i-1][j][k];
			double fy_jph = (j==jsize-1) ? 0 : f[i][j+1][k]-f000, fy_jmh = (j==0) ? 0 : f000-f[i][j-1][k];
			double fz_kph = (i==isize-1) ? 0 : f[i][j][k+1]-f000, fz_kmh = (k==0) ? 0 : f000-f[i][j][k-1];
			fxx[i][j][k] = fx_iph-fx_imh;
			fyy[i][j][k] = fy_jph-fy_jmh; 
			fzz[i][j][k] = fz_kph-fz_kmh;}}
	// interpolation
	double operator()( double x, double y, double z ) const {
		// physical coordinates => logical coordinates
		double id = grid.i_fr_x(x); int i = (int)id; if(i<0)i=0; if(i+1>=grid.isize) i=grid.isize-2;
		double jd = grid.j_fr_y(y); int j = (int)jd; if(j<0)j=0; if(j+1>=grid.jsize) j=grid.jsize-2;
		double kd = grid.k_fr_z(z); int k = (int)kd; if(k<0)k=0; if(k+1>=grid.ksize) k=grid.ksize-2;
		id-=i;
		jd-=j;
		kd-=k;
		const ARRAY3D<double>& f = *pf;
		double trilinear =(1-id)*((1-jd)*((1-kd)*f[i  ][j  ][k  ]+(kd)*f[i  ][j  ][k+1])
		                         +(  jd)*((1-kd)*f[i  ][j+1][k  ]+(kd)*f[i  ][j+1][k+1]))
		                 +(  id)*((1-jd)*((1-kd)*f[i+1][j  ][k  ]+(kd)*f[i+1][j  ][k+1])
				                 +(  jd)*((1-kd)*f[i+1][j+1][k  ]+(kd)*f[i+1][j+1][k+1]));
		// quadratic updates
		double fxx_h00 = MINMOD(fxx[i][j  ][k  ],fxx[i+1][j  ][k  ]);
		double fxx_h01 = MINMOD(fxx[i][j  ][k+1],fxx[i+1][j  ][k+1]);
		double fxx_h10 = MINMOD(fxx[i][j+1][k  ],fxx[i+1][j+1][k  ]);
		double fxx_h11 = MINMOD(fxx[i][j+1][k+1],fxx[i+1][j+1][k+1]);
		double fyy_0h0 = MINMOD(fyy[i  ][j][k  ],fyy[i  ][j+1][k  ]);
		double fyy_0h1 = MINMOD(fyy[i  ][j][k+1],fyy[i  ][j+1][k+1]);
		double fyy_1h0 = MINMOD(fyy[i+1][j][k  ],fyy[i+1][j+1][k  ]);
		double fyy_1h1 = MINMOD(fyy[i+1][j][k+1],fyy[i+1][j+1][k+1]);
		double fzz_00h = MINMOD(fzz[i  ][j  ][k],fzz[i  ][j  ][k+1]);
		double fzz_01h = MINMOD(fzz[i  ][j+1][k],fzz[i  ][j+1][k+1]);
		double fzz_10h = MINMOD(fzz[i+1][j  ][k],fzz[i+1][j  ][k+1]);
		double fzz_11h = MINMOD(fzz[i+1][j+1][k],fzz[i+1][j+1][k+1]);
		double fxx = (1-jd)*((1-kd)*fxx_h00+(kd)*fxx_h01)
				   + (  jd)*((1-kd)*fxx_h10+(kd)*fxx_h11);
		double fyy = (1-id)*((1-kd)*fyy_0h0+(kd)*fyy_0h1)
				   + (  id)*((1-kd)*fyy_1h0+(kd)*fyy_1h1);
		double fzz = (1-id)*((1-jd)*fzz_00h+(jd)*fzz_01h)
				   + (  id)*((1-jd)*fzz_10h+(jd)*fzz_11h);
		return trilinear- .5*fxx*id*(1-id)
			            - .5*fyy*jd*(1-jd)
						- .5*fzz*kd*(1-kd);
	}
};

#endif
