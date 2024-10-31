#ifndef QUADRATIC_2D_H
#define QUADRATIC_2D_H

#include "../Macros.h"
#include "GRID2D.h"
#include "../ARRAY2D.h"

class QuadraticInterpolation2D : public Function_2D {
protected:
	GRID2D grid;
	const ARRAY2D<double>* pf;
	      ARRAY2D<double> fxx,fyy;
public:
	QuadraticInterpolation2D( const GRID2D& grid_, const ARRAY2D<double>& f ) {
		grid=grid_; int isize=grid.isize;
		pf=&f;      int jsize=grid.jsize;
		// 
		fxx.resize(isize,jsize);
		fyy.resize(isize,jsize);
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++) {
			double fx_iph = (i==isize-1) ? 0 : f[i+1][j]-f[i][j], fx_imh = (i==0) ? 0 : f[i][j]-f[i-1][j];
			double fy_jph = (j==jsize-1) ? 0 : f[i][j+1]-f[i][j], fy_jmh = (j==0) ? 0 : f[i][j]-f[i][j-1];
			fxx[i][j] = fx_iph-fx_imh;
			fyy[i][j] = fy_jph-fy_jmh; }
	}
	// interpolation
	double operator()( double x, double y ) const {
		// physical coordinates => logical coordinates
		double id = grid.i_fr_x(x); int i = (int)id; if(i<0)i=0; if(i+1>=grid.isize) i=grid.isize-2;
		double jd = grid.j_fr_y(y); int j = (int)jd; if(j<0)j=0; if(j+1>=grid.jsize) j=grid.jsize-2;
		id-=i;
		jd-=j;
		const ARRAY2D<double>& f = *pf;
		double bilinear = (1-id)*((1-jd)*f[i  ][j]+(jd)*f[i  ][j+1])
		                 +(  id)*((1-jd)*f[i+1][j]+(jd)*f[i+1][j+1]);
		// quadratic updates
		double fxx_h0 = MINMOD(fxx[i][j  ],fxx[i+1][j  ]);
		double fxx_h1 = MINMOD(fxx[i][j+1],fxx[i+1][j+1]);
		double fyy_0h = MINMOD(fyy[i  ][j],fyy[i  ][j+1]);
		double fyy_1h = MINMOD(fyy[i+1][j],fyy[i+1][j+1]);
		double fxx = fxx_h0*(1-jd) + fxx_h1*(jd);
		double fyy = fyy_0h*(1-id) + fyy_1h*(id);
		return bilinear - .5*fxx*id*(1-id)
			            - .5*fyy*jd*(1-jd);
	}
};

#endif
